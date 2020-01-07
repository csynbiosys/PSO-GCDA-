from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField
from wtforms.validators import DataRequired

class SpecificationForm(FlaskForm):
    circuit_name = StringField('Circuit Name: ', validators=[DataRequired()])
    gene_nodes = StringField('Genes:', validators=[DataRequired()])
    input_nodes = StringField('Inputs:', validators=[DataRequired()])
    output_nodes = StringField('Outputs:', validators=[DataRequired()])
    max_expression = StringField('Maximum Protein Expression:', validators=[DataRequired()])
    time_span = StringField('Time Span:', validators=[DataRequired()])
    submit = SubmitField('Submit')
