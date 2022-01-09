from math import sqrt
import matplotlib.pyplot as plt
import json
import skfuzzy as fuzz
from skfuzzy import control as ctrl
import numpy as np

class UAR:
    pid = {
        'error': [0, ],             #uchyb regulacji
        'gain': 1.0,                #wartość wzmocnienia regulatora
        'sample_time': 0.05,        #czas próbkowania
        'differential_time': 0.05,  #czas wyprzedzenia
        'integration_time': 0.75,   #czas zdwojenia
        'h_z': [],                  #wartość zadana
        'h': [0, ]                  #poziom substancji w zbiorniku
    }
    valvee = {
        'u_max': 10,                #wartość max poziomu substancji w zbiorniku
        'u_min': -10,               #wartość min poziomu substancji w zbiorniku
        'u': [0.1, ],               #wartość aktualna sygnału sterującego
        'Q_d_min': 0,               #minimalne natężenie dopływu
        'Q_d_max': 1,               #maksymalne natężenie dopływu
        'Q_d': []                   #wartość aktualna natężenia dopływu
    }
    tank = {
        'h_max': 10,                #maksymalny poziom substancji w zbiorniku
        'h_min': 0,                 #minimalny poziom substancji w zbiorniku
        'A': 2.5,                   #pole powierzchni przekroju poprzecznego zbiornika
        'B': 0.25,                  #współczynnik wypływu
        'Q_o': [],                  #natężenie odpływu
    }
    N = 10000
    x = [0, ]


    def JsonDate(self):
        listatoUaR = []
        with open('data.json') as json_data:
            data_dict = json.load(json_data)
            for key, value in data_dict.items():
                listatoUaR.append(value)
        listatoUaR.remove(listatoUaR[0])
        listatoUaR.remove(listatoUaR[-1])

        self.pid['sample_time'] = float(listatoUaR[0])
        self.pid['differential_time'] = float(listatoUaR[1])
        self.pid['integration_time'] = float(listatoUaR[2])
        self.pid['gain'] = float(listatoUaR[3])
        self.pid['h_z'].append(float(listatoUaR[4]))
        self.tank['A'] = float(listatoUaR[5])
        self.tank['B'] = float(listatoUaR[6])
        self.tank['h_max'] = int(listatoUaR[7])
        self.valvee['u_max'] = int(listatoUaR[8])
        self.valvee['u_min'] = int(listatoUaR[9])
        self.valvee['Q_d_max'] = int(listatoUaR[10])

    def reset_data(self):

        self.pid['h_z'] = []
        self.pid['h'] = [0, ]
        self.N = 10000
        self.x = [0, ]
        plt.clf()

    # PID
    def pid_controler(self):
        self.pid['error'].append((self.pid['h_z'][-1] - self.pid['h'][-1]) / 10)
        self.valvee['u'].append(
            (self.pid['gain'] * (self.pid['error'][-1] + (self.pid['sample_time'] / self.pid['integration_time']) * sum(
                self.pid['error']) + (self.pid['differential_time'] / self.pid['sample_time']) * (
                                             self.pid['error'][1] - self.pid['error'][-1]))) / 10)

    # Zawór
    def valve(self):
        if (self.valvee['u'][-1] <= self.valvee['Q_d_min']):
            self.valvee['Q_d'].append(self.valvee['u_min'])
        elif (self.valvee['u'][-1] >= self.valvee['Q_d_max']):
            self.valvee['Q_d'].append(self.valvee['u_max'])
        else:
            self.valvee['Q_d'].append(self.valvee['u'][-1])

    # Zbiornik
    def tankk(self):
        self.tank['Q_o'].append(self.tank['B'] * sqrt(self.pid['h'][-1]))
        self.pid['h'].append(max(min(((-self.tank['Q_o'][-1] + self.valvee['Q_d'][-1]) * self.pid['sample_time'] /
                                      self.tank['A'] + self.pid['h'][-1]), self.tank['h_max']), self.tank['h_min']))

    def count(self):
        self.reset_data()
        self.JsonDate()
        for n in range(0, self.N):
            self.x.append(obiekt.pid['sample_time'] * n)
            self.pid['h_z'].append(obiekt.pid['h_z'][-1])
            self.pid_controler()
            self.valve()
            self.tankk()

    def display(self):
        self.count()
        plt.xlabel('t [s]')
        plt.ylabel('h [m]')
        plt.title('Wykres zależności poziomu substancji w zbiorniku', fontsize=14)
        plt.plot(self.x, self.pid['h'], label='$h$', color='blue')
        plt.plot(self.x, self.pid['h_z'], label='$h_z$', color='red', linewidth='2', linestyle='--')
        plt.legend()
        plt.savefig('static/plot.png')
# #FUZZY
#     def Mlr(self):
#         x_e = ctrl.Antecedent(np.arange(-1, 1, 0.00001), 'x_e')
#         x_ce = ctrl.Antecedent(np.arange(-1, 1, 0.00001), 'x_ce')
#         y_du = ctrl.Consequent(np.arange(-1, 1, 0.00001), 'y_du')
#
#         # x_e
#         x_e['e_DU'] = fuzz.trimf(x_e.universe, [-1.3, -1, -0.667])
#         x_e['e_SU'] = fuzz.trimf(x_e.universe, [-1, -0.667, -0.333])
#         x_e['e_MU'] = fuzz.trimf(x_e.universe, [-0.667, -0.333, 0.0])
#         x_e['e_Z'] = fuzz.trimf(x_e.universe, [-0.333, 0.0, 0.333])
#         x_e['e_MD'] = fuzz.trimf(x_e.universe, [0.0, 0.333, 0.667])
#         x_e['e_SD'] = fuzz.trimf(x_e.universe, [0.333, 0.667, 1])
#         x_e['e_DD'] = fuzz.trimf(x_e.universe, [0.667, 1, 1.33])
#
#         # x_ce
#         x_ce['ce_DU'] = fuzz.trimf(x_ce.universe, [-1.3, -1, -0.667])
#         x_ce['ce_SU'] = fuzz.trimf(x_ce.universe, [-1, -0.667, -0.333])
#         x_ce['ce_MU'] = fuzz.trimf(x_ce.universe, [-0.667, -0.333, 0.0])
#         x_ce['ce_Z'] = fuzz.trimf(x_ce.universe, [-0.333, 0.0, 0.333])
#         x_ce['ce_MD'] = fuzz.trimf(x_ce.universe, [0.0, 0.333, 0.667])
#         x_ce['ce_SD'] = fuzz.trimf(x_ce.universe, [0.333, 0.667, 1])
#         x_ce['ce_DD'] = fuzz.trimf(x_ce.universe, [0.667, 1, 1.33])
#
#         # y_du
#         y_du['du_BDU'] = fuzz.trimf(y_du.universe, [-1.25, -1, -0.75])
#         y_du['du_DU'] = fuzz.trimf(y_du.universe, [-1, -0.75, -0.5])
#         y_du['du_SU'] = fuzz.trimf(y_du.universe, [-0.75, -0.5, -0.25])
#         y_du['du_MU'] = fuzz.trimf(y_du.universe, [-0.5, -0.25, 0])
#         y_du['du_Z'] = fuzz.trimf(y_du.universe, [-0.25, 0, 0.25])
#         y_du['du_MD'] = fuzz.trimf(y_du.universe, [0, 0.25, 0.5])
#         y_du['du_SD'] = fuzz.trimf(y_du.universe, [0.25, 0.5, 0.75])
#         y_du['du_DD'] = fuzz.trimf(y_du.universe, [0.5, 0.75, 1])
#         y_du['du_BDD'] = fuzz.trimf(y_du.universe, [0.75, 1, 1.25])
#
#         rule1 = ctrl.Rule(x_e['e_DU'] & x_ce['ce_DU'], y_du['du_BDU'])
#         rule2 = ctrl.Rule(x_e['e_DU'] & x_ce['ce_SU'], y_du['du_BDU'])
#         rule3 = ctrl.Rule(x_e['e_DU'] & x_ce['ce_MU'], y_du['du_BDU'])
#         rule4 = ctrl.Rule(x_e['e_DU'] & x_ce['ce_Z'], y_du['du_DU'])
#         rule5 = ctrl.Rule(x_e['e_DU'] & x_ce['ce_MD'], y_du['du_SU'])
#         rule6 = ctrl.Rule(x_e['e_DU'] & x_ce['ce_SD'], y_du['du_MU'])
#         rule7 = ctrl.Rule(x_e['e_DU'] & x_ce['ce_DD'], y_du['du_Z'])
#
#         rule8 = ctrl.Rule(x_e['e_SU'] & x_ce['ce_DU'], y_du['du_BDU'])
#         rule9 = ctrl.Rule(x_e['e_SU'] & x_ce['ce_SU'], y_du['du_BDU'])
#         rule10 = ctrl.Rule(x_e['e_SU'] & x_ce['ce_MU'], y_du['du_DU'])
#         rule11 = ctrl.Rule(x_e['e_SU'] & x_ce['ce_Z'], y_du['du_SU'])
#         rule12 = ctrl.Rule(x_e['e_SU'] & x_ce['ce_MD'], y_du['du_MU'])
#         rule13 = ctrl.Rule(x_e['e_SU'] & x_ce['ce_SD'], y_du['du_Z'])
#         rule14 = ctrl.Rule(x_e['e_SU'] & x_ce['ce_DD'], y_du['du_MD'])
#
#         rule15 = ctrl.Rule(x_e['e_MU'] & x_ce['ce_DU'], y_du['du_BDU'])
#         rule16 = ctrl.Rule(x_e['e_MU'] & x_ce['ce_SU'], y_du['du_DU'])
#         rule17 = ctrl.Rule(x_e['e_MU'] & x_ce['ce_MU'], y_du['du_SU'])
#         rule18 = ctrl.Rule(x_e['e_MU'] & x_ce['ce_Z'], y_du['du_MU'])
#         rule19 = ctrl.Rule(x_e['e_MU'] & x_ce['ce_MD'], y_du['du_Z'])
#         rule20 = ctrl.Rule(x_e['e_MU'] & x_ce['ce_SD'], y_du['du_MD'])
#         rule21 = ctrl.Rule(x_e['e_MU'] & x_ce['ce_DD'], y_du['du_SD'])
#
#         rule22 = ctrl.Rule(x_e['e_Z'] & x_ce['ce_DU'], y_du['du_DU'])
#         rule23 = ctrl.Rule(x_e['e_Z'] & x_ce['ce_SU'], y_du['du_SU'])
#         rule24 = ctrl.Rule(x_e['e_Z'] & x_ce['ce_MU'], y_du['du_MU'])
#         rule25 = ctrl.Rule(x_e['e_Z'] & x_ce['ce_Z'], y_du['du_Z'])
#         rule26 = ctrl.Rule(x_e['e_Z'] & x_ce['ce_MD'], y_du['du_MD'])
#         rule27 = ctrl.Rule(x_e['e_Z'] & x_ce['ce_SD'], y_du['du_SD'])
#         rule28 = ctrl.Rule(x_e['e_Z'] & x_ce['ce_DD'], y_du['du_DD'])
#
#         rule29 = ctrl.Rule(x_e['e_MD'] & x_ce['ce_DU'], y_du['du_SU'])
#         rule30 = ctrl.Rule(x_e['e_MD'] & x_ce['ce_SU'], y_du['du_MU'])
#         rule31 = ctrl.Rule(x_e['e_MD'] & x_ce['ce_MU'], y_du['du_Z'])
#         rule32 = ctrl.Rule(x_e['e_MD'] & x_ce['ce_Z'], y_du['du_MD'])
#         rule33 = ctrl.Rule(x_e['e_MD'] & x_ce['ce_MD'], y_du['du_SD'])
#         rule34 = ctrl.Rule(x_e['e_MD'] & x_ce['ce_SD'], y_du['du_DD'])
#         rule35 = ctrl.Rule(x_e['e_MD'] & x_ce['ce_DD'], y_du['du_BDD'])
#
#         rule36 = ctrl.Rule(x_e['e_SD'] & x_ce['ce_DU'], y_du['du_MU'])
#         rule37 = ctrl.Rule(x_e['e_SD'] & x_ce['ce_SU'], y_du['du_Z'])
#         rule38 = ctrl.Rule(x_e['e_SD'] & x_ce['ce_MU'], y_du['du_MD'])
#         rule39 = ctrl.Rule(x_e['e_SD'] & x_ce['ce_Z'], y_du['du_SD'])
#         rule40 = ctrl.Rule(x_e['e_SD'] & x_ce['ce_MD'], y_du['du_DD'])
#         rule41 = ctrl.Rule(x_e['e_SD'] & x_ce['ce_SD'], y_du['du_BDD'])
#         rule42 = ctrl.Rule(x_e['e_SD'] & x_ce['ce_DD'], y_du['du_BDD'])
#
#         rule43 = ctrl.Rule(x_e['e_DD'] & x_ce['ce_DU'], y_du['du_Z'])
#         rule44 = ctrl.Rule(x_e['e_DD'] & x_ce['ce_SU'], y_du['du_MD'])
#         rule45 = ctrl.Rule(x_e['e_DD'] & x_ce['ce_MU'], y_du['du_SD'])
#         rule46 = ctrl.Rule(x_e['e_DD'] & x_ce['ce_Z'], y_du['du_DD'])
#         rule47 = ctrl.Rule(x_e['e_DD'] & x_ce['ce_MD'], y_du['du_BDD'])
#         rule48 = ctrl.Rule(x_e['e_DD'] & x_ce['ce_SD'], y_du['du_BDD'])
#         rule49 = ctrl.Rule(x_e['e_DD'] & x_ce['ce_DD'], y_du['du_BDD'])
#
#         global fuzzy_ctrl
#         fuzzy_ctrl = ctrl.ControlSystem(
#             [rule1, rule2, rule3, rule4, rule5, rule6, rule7, rule8, rule9, rule10, rule11, rule12,
#              rule13, rule14, rule15, rule16, rule17, rule18, rule19, rule20, rule21, rule22, rule23,
#              rule24, rule25, rule26, rule27, rule28, rule29, rule30, rule31, rule32, rule33, rule34,
#              rule35, rule36, rule37, rule38, rule39, rule40, rule41, rule42, rule43, rule44, rule45,
#              rule46, rule47, rule48, rule49])
#
#         fuzzy_simulation = ctrl.ControlSystemSimulation(fuzzy_ctrl)
#
#     def fuzzy(u, V_min, V_max, nr):
#
#         k_e = 12.0
#         k_c = 12.0
#         u1 = k_e * u
#         u2 = (u - u_last_fuzzy[nr]) * 1 / k_c
#         u_last_fuzzy[nr] = u
#
#         if u1 < -1.0:
#             u1 = -1.0
#         if u1 > 1.0:
#             u1 = 1.0
#
#         if u2 < -1.0:
#             u2 = -1.0
#         if u2 > 1.0:
#             u2 = 1.0
#
#         fuzzy_simulation.input['x_e'] = u1
#         fuzzy_simulation.input['x_ce'] = u2
#         fuzzy_simulation.compute()
#
#         print(fuzzy_simulation.output['y_du'])
#         # y_du.view(sim=fuzzy_simulation)
#
#         yreg_fuzzy[nr] = fuzzy_simulation.output['y_du']
#
#         if yreg_fuzzy[nr] < V_min:
#             yreg_fuzzy[nr] = V_min
#         if yreg_fuzzy[nr] > V_max:
#             yreg_fuzzy[nr] = V_max
#
#         return yreg_fuzzy[nr] * 2.0
obiekt = UAR()