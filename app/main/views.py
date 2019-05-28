from flask import render_template, session, redirect, url_for, current_app, flash
# from .. import db
# from ..models import User
# from ..email import send_email
from . import main
from .forms import DataForm

import models

@main.route('/', methods=['GET', 'POST'])
def index():
    DataForm.init(session.get('sequences'))
    form = DataForm()
    if form.validate_on_submit():
        # user = User.query.filter_by(username=form.name.data).first()
        # if user is None:
        #     user = User(username=form.name.data)
        #     db.session.add(user)
        #     db.session.commit()
        #     session['known'] = False
        #     if current_app.config['FLASKY_ADMIN']:
        #         send_email(current_app.config['FLASKY_ADMIN'], 'New User',
        #                    'mail/new_user', user=user)
        # else:
        #     session['known'] = True
        session['sequences'] = form.sequences.data
        # print(session['sequences'])
        # print(type(session['sequences']))
        # input data is a string separated by \n

        # Do some processing here.
        X, processing_error = models.preprocess_seq(form.sequences.data)

        # debug
        if X is not None:
            print('X shape = ', X.shape)
        print('Error = ', processing_error)

        if processing_error:
           flash(processing_error)
           session['results'] = []
        else:
            # Save results in session.
            session['results'] = models.identify_enhancer(X)
            print(session['results'])
            # session['results'] = [('enhancer', 'strong'), ('enhancer', 'weak'), ('non-enhancer', 'non-enhancer'), ('non-enhancer', 'non-enhancer'),
            #                       ('enhancer', 'strong'), ('enhancer', 'weak'), ('non-enhancer', 'non-enhancer')]

        return redirect(url_for('.index'))
    return render_template('index.html', form=form, results=session.get('results'))