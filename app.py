from flask import Flask, render_template, request
from flask_wtf import FlaskForm
from wtforms import FloatField, DecimalRangeField, SubmitField
from wtforms.validators import InputRequired
import json
from compute import UAR, obiekt


app = Flask(__name__)
app.config['SECRET_KEY'] = 'mykey'
app.config['DEBUG'] = True

class MODEL(FlaskForm):
    default_value = [0.05, 0.05, 0.75, 1, 3, 2.5, 0.25, 10, 10, -10, 1, 0.20, 0.03]
    sample_time = DecimalRangeField('Okres próbkowania: ', default=default_value[0])
    differential_time = DecimalRangeField('Czas wyprzedzenia: ', default=default_value[1])
    integration_time = DecimalRangeField('Czas zdwojenia :', default=default_value[2])
    gain = DecimalRangeField('Wartość wzmocnienia regulatora: ', default=default_value[3])
    h_z = DecimalRangeField('Wartość zadana: ', default=default_value[4])
    A = FloatField('Pole powierzchni przekroju poprzecznego: ', default=default_value[5], validators=[InputRequired()])
    B = FloatField('Współczynnik wypływu: ', default=default_value[6], validators=[InputRequired()])
    h_max = FloatField('Maksymalny poziom substancji w zbiorniku: ', default=default_value[7], validators=[InputRequired()])
    u_max = FloatField('Maksymalna wartość sygnału sterującego: ', default=default_value[8], validators=[InputRequired()])
    u_min = FloatField('Minimalna wartość sygnału sterującego: ', default=default_value[9], validators=[InputRequired()])
    Q_d_max = FloatField('Maksymalne natężenie dopływu: ', default=default_value[10], validators=[InputRequired()])
    error = FloatField('Fuzzy Error', default=default_value[11], validators=[InputRequired()])
    errorChange = FloatField('Fuzzy Error Change', default=default_value[12], validators=[InputRequired()])
    submit_all = SubmitField('Zatwierdź Dane')

@app.route('/', methods = ["POST","GET"])
def index():
    form = MODEL(request.form)

    if request.method == "POST" and form.validate():
        with open('data.json', 'w') as f:
            json.dump(request.form, f)
    return render_template('view_pid_simulator.html', form=form, home_active="class=active" )

@app.route("/Generate_plot")
def SomeFunction():
    UAR.display(obiekt)
    return "plot.png"

@app.route("/Generate_plot_Fuzz")
def GenerateFuzzy():
    UAR.FuzzyDisplay(obiekt)
    return "Fuzze.png"

@app.route("/Generate_porownanie")
def GenerateSet():
    UAR.Porownanie(obiekt)
    return "porownanie.png"

@app.route("/AboutMe")
def GetAboutMe():
    return render_template('AboutMe.html', title="O Autorze", AboutMe_active="class=active")

@app.route("/AboutPID")
def GetAboutPID():
    return render_template('AboutPID.html', title="PID", PID_active="class=active")

@app.route("/AboutFuzzyPID")
def GetAboutFuzzyPID():
    return render_template('AboutFuzzyPID.html', title="Fuzzy PID", FuzzyPID_active="class=active")


if __name__ == '__main__':
    app.run()
