from flask import Flask, render_template, flash, request
from forms import LoginForm
from config import Config
import logging
import vcf
import numpy as np
from utils import *
import folium

app = Flask(__name__)
app.config.from_object(Config)


@app.route('/', methods=['GET', 'POST'])
def index():
    form = LoginForm()
    # coordinates for Republic of Belarus
    start_coords = (53.9211011, 27.568865)
    folium_map = folium.Map(location=start_coords, zoom_start=6, tiles='cartodbpositron')
    request_haplogroup = None

    if request.method == "POST":

        if request.files:
            image = request.files["file"]
            if allowed_image(image.filename, app.config["ALLOWED_IMAGE_EXTENSIONS"]):
                fb = image.read()
                vcf_reader = vcf.Reader(fb.decode().split("\n"))
                all_groups, group = get_groups_from_vcf(app.config["DEPENDENCIES_FOLDER"], vcf_reader, mapped_markers,
                                                        y_tree)
                print(group)
                form.haplogroup.data = group
                form.tree.label = " -> ".join(all_groups)
                request_haplogroup = group
            else:
                logging.info('Sorry, only vcf format is supported')

    if form.validate_on_submit():
        request_haplogroup = form.haplogroup.data.strip()

    if request_haplogroup is not None:
        if request_haplogroup in frequency_table.columns:
            df = get_probability_of_population(request_haplogroup, frequency_table)

            # Add the color for the chloropleth:
            folium_map.choropleth(
                geo_data=gj1,  # locations_df,
                name='choropleth',
                data=df,  # state_data,
                columns=['id_region', request_haplogroup],
                key_on='feature.id',
                fill_color="GnBu",  # "YlGn"
                fill_opacity=0.7,
                line_opacity=0.8,
                nan_fill_color='#ffffff00'
            )
            folium.LayerControl().add_to(folium_map)
            logging.info('The requested haplogroup is {}'.format(request_haplogroup))
        else:
            logging.info('Sorry, there is no data about {} haplogroup'.format(request_haplogroup))
            form.exception.label = 'Sorry, there is no data about {} haplogroup'.format(request_haplogroup)

    folium_map.save('templates/map.html')
    return render_template('index.html', form=form)


global frequency_table, gj1, mapped_markers, y_tree


def transform_frequency_table(df, reg_dict):
    df = df.T
    df["rus_region"] = df.index
    df["id_region"] = df["rus_region"].apply(
        lambda x: reg_dict[x] if x in reg_dict else np.NaN)
    df.fillna(0, inplace=True)
    return df


if __name__ == '__main__':
    logging.basicConfig(filename="logfile.log", level=logging.INFO)
    existed_groups = set()
    probabilty_dict = {"total_amount": {}}

    logging.info("Program started upload dependent files")
    try:
        mapped_markers, y_tree, gj1, regions_dict = load_dependencies(app.config["DEPENDENCIES_FOLDER"])
    except:
        logging.exception("Can't read dependencies. Check directories, please")
        raise

    logging.info("Dependent files were sucсessfully downloaded")
    logging.info("Program started upload data")
    for file in os.listdir(app.config["DATA_FOLDER"]):
        excel_file_pth = os.path.join(app.config["DATA_FOLDER"], file)
        batch_data = pd.read_excel(excel_file_pth)
        chrY_data = filter_column_by_value(batch_data, "#CHROM", "chrY")
        drop_females(chrY_data)
        group_haplo_dict = get_group_and_branch(chrY_data, y_tree, mapped_markers)
        add_existed_group(group_haplo_dict, existed_groups)
        location_df = get_location_mapping(excel_file_pth)
        probabilty_dict = create_population_frequency_table(probabilty_dict, group_haplo_dict,
                                                            location_df, existed_groups)
    logging.info("The data was downloaded")

    frequency_table = pd.DataFrame.from_dict(probabilty_dict, orient='index')
    frequency_table = transform_frequency_table(frequency_table, regions_dict)
    app.run(debug=True)
