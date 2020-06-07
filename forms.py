from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, TextField
from wtforms.validators import DataRequired


class LoginForm(FlaskForm):
    haplogroup = StringField('Input haplogroup', validators=[DataRequired()])
    tree = TextField("")
    submit = SubmitField('Find regions')
    exception = StringField("")

