import os
import json
import geojson
import pandas as pd


def get_markers_from_vcf_file(pth, vcf_reader, marker_mapping):
    df = pd.read_excel(os.path.join(pth, "marker_mapping.xlsx"))
    existed_pos = set(df.loc[:, "based on GRCh38/hg38"])
    markers = []

    for record in vcf_reader:
        if record.POS in existed_pos:
                marker = df[df.loc[:, "based on GRCh38/hg38"] == record.POS]["Marker name"].values[0]
                base = record.ALT[0] if record.ALT[0] is not None else record.REF
                if marker in marker_mapping and marker_mapping[marker] == base:
                    markers.append(marker)
    return markers


def get_groups_from_vcf(pth, vcf_reader, marker_mapping, y_tree):
    markers = get_markers_from_vcf_file(pth, vcf_reader, marker_mapping)
    all_groups, group = get_haplogroup_by_markers(y_tree, set(markers), groups=[], previous_group=None)
    return all_groups, group
    
    
def allowed_image(filename, allowed_extensions):

    # We only want files with a . in the filename
    file_extension = os.path.splitext(filename)[1]

    # Check if the extension is in ALLOWED_IMAGE_EXTENSIONS
    if file_extension.upper() in allowed_extensions:
        return True
    else:
        return False
        
        
def get_location_mapping(excel_file_pth):
    df = pd.read_excel(excel_file_pth, sheet_name="mapping")
    df.dropna(inplace=True)
    df.set_index("Sample", inplace=True)
    return df


def load_dependencies(pth):
    with open(os.path.join(pth, "markers_map.json"), "r") as f:
        mapped_markers = json.load(f)

    with open(os.path.join(pth, "y_tree.json"), "r") as f:
        y_tree = json.load(f)

    path_to_geo_file = os.path.join(pth, "geo.geojson")
    with open(path_to_geo_file) as f:
        gj1 = geojson.load(f)

    with open(os.path.join(pth, "regions_dictionary.json"), "r", encoding="utf-8") as f:
        regions_dict = json.load(f)
    return mapped_markers, y_tree, gj1, regions_dict


def read_batch_files_to_dataframe(dir_to_batch_files):
    batch_data = pd.DataFrame()
    for file in os.listdir(dir_to_batch_files):
        batch_data = batch_data.append(pd.read_excel(os.path.join(dir_to_batch_files, file)))
    return batch_data


def filter_column_by_value(df, column, value):
    return df[df[column] == value]


def drop_females(df):
    for column in df.columns:
        if all(df[column].isna()):
            df.drop(columns=[column], inplace=True)


def get_markers_from_column(df, sample_name_column, marker_name_column, mapped_markers):
    sample_markers = []
    for ind in df.index:
        marker = df.loc[ind, marker_name_column]
        if marker in mapped_markers:
            if mapped_markers[marker] == df[sample_name_column][ind]:
                sample_markers.append(marker)
    return sample_markers


def get_haplogroup_by_markers(y_tree, set_sample_markers, groups=[], previous_group=None):
    children = {}
    for group in y_tree:
        common_markers = set(y_tree[group]['markers']) & set_sample_markers

        if len(common_markers) != 0:
            for shared_mrk in common_markers:
                set_sample_markers.remove(shared_mrk)

            if len(y_tree[group]['descendants']) == 0:
                groups.append(group)
                return groups, group

            else:
                if "parent" in y_tree[group]:
                    groups += y_tree[group]['parent']
                groups.append(group)
                return get_haplogroup_by_markers(y_tree[group]['descendants'], set_sample_markers, groups, group)

        elif len(y_tree[group]['descendants']) != 0:
            for descendant in y_tree[group]['descendants']:
                new_descendant = y_tree[group]['descendants'][descendant].copy()
                if "parent" in y_tree[group]:
                    new_descendant["parent"] = y_tree[group]["parent"].copy()
                    new_descendant["parent"].append(group)
                else:
                    new_descendant["parent"] = [group]
                children[descendant] = new_descendant

    if len(children):
        return get_haplogroup_by_markers(children, set_sample_markers, groups, previous_group)
    else:
        return groups, previous_group


def get_group_and_branch(chrY_data, y_tree, mapped_markers):
    group_haplo_dict = {}
    for sample_name_column in chrY_data.columns[9:]:
        sample_markers_list = get_markers_from_column(chrY_data, sample_name_column, marker_name_column="ANNOTATION",
                                                      mapped_markers=mapped_markers)
        list_of_branches, haplogroup = get_haplogroup_by_markers(y_tree.copy(), set(sample_markers_list), [])
        print(sample_name_column, list_of_branches, haplogroup)
        group_haplo_dict[sample_name_column] = {"group": haplogroup, "all_groups": list_of_branches}
        # print(list_of_branches)
    return group_haplo_dict


def add_existed_group(group_haplo_dict, existed_groups):
    for sample in group_haplo_dict:
        for group in group_haplo_dict[sample]["all_groups"]:
            existed_groups.add(group)


def increment_dict_value(dictionary, key, value=1):
    if key in dictionary:
        dictionary[key] += value
    else:
        dictionary[key] = value


'''def create_population_frequency_table_1(group_haplo_dict):
    probabilty_dict = {"total_amount": {}}
    for sample in group_haplo_dict:

        if sample in location_df.index:
            # print(sample)
            sample_population = location_df.loc[sample]["population_id"]
            increment_dict_value(probabilty_dict["total_amount"], sample_population)
            # probabilty_dict["total_amount"][sample_population] += 1
            for marker in group_haplo_dict[sample]["all_groups"]:
                if marker not in probabilty_dict:
                    probabilty_dict[marker] = {}
                # print(marker)
                increment_dict_value(probabilty_dict[marker], sample_population)

    return probabilty_dict
'''


def create_population_frequency_table(probabilty_dict, group_haplo_dict, location_df, existed_groups):
    for sample in group_haplo_dict:

        if sample in location_df.index:
            # print(sample)
            sample_population = location_df.loc[sample]["population_id"]
            increment_dict_value(probabilty_dict["total_amount"], sample_population)
            # probabilty_dict["total_amount"][sample_population] += 1
            for group in existed_groups:
                if group not in probabilty_dict:
                    probabilty_dict[group] = {}
                if group in group_haplo_dict[sample]["all_groups"]:
                    # print(marker)
                    increment_dict_value(probabilty_dict[group], sample_population)
                else:
                    increment_dict_value(probabilty_dict[group], sample_population, value=0)
    return probabilty_dict


def get_probability_of_population(haplogroup, frequency_table):
    df = frequency_table[[haplogroup, "total_amount", "rus_region", "id_region"]]
    df[haplogroup] /= df["total_amount"]
    #df["region"] = df["index"].apply(lambda x: reg_dict[x] if x in reg_dict else np.NaN)
    return df