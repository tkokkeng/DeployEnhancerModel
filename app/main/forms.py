from flask_wtf import FlaskForm
from wtforms import TextAreaField, SubmitField
from wtforms.validators import DataRequired

class DataForm(FlaskForm):
    # sequences = TextAreaField('Sequence Data', validators=[DataRequired()])
    # submit = SubmitField('Submit')
    # sequences = None
    # submit = None

    @classmethod
    def init(cls, data):
        DataForm.sequences = TextAreaField('Sequence Data', validators=[DataRequired()], default=data,
                                           render_kw={"rows": 10, "cols": 120})
        DataForm.submit = SubmitField('Submit')
