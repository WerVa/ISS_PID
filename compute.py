from math import sqrt
import matplotlib.pyplot as plt
import json
import skfuzzy as fuzz
from skfuzzy import control as ctrl
import numpy as np

class UAR:
    pid = {
        'error': [0, ],             # uchyb regulacji
        'gain': 1.0,                # wartość wzmocnienia regulatora
        'sample_time': 0.05,        # czas próbkowania
        'differential_time': 0.05,  # czas wyprzedzenia
        'integration_time': 0.75,   # czas zdwojenia
        'h_z': [],                  # wartość zadana
        'h': [0, ]                  # poziom substancji w zbiorniku
    }
    valvee = {
        'u_max': 10,                # wartość max poziomu substancji w zbiorniku
        'u_min': -10,               # wartość min poziomu substancji w zbiorniku
        'u': [0.1, ],               # wartość aktualna sygnału sterującego
        'Q_d_min': 0,               # minimalne natężenie dopływu
        'Q_d_max': 1,               # maksymalne natężenie dopływu
        'Q_d': [],                  # wartość aktualna natężenia dopływu
        'h_start': 0,               # wartość początkowa
        'y': [[], [], [], [], []],  # macierz dla h, Q_d, Q_o, h_z

    }
    tank = {
        'h_max': 10,                # maksymalny poziom substancji w zbiorniku
        'h_min': 0,                 # minimalny poziom substancji w zbiorniku
        'A': 2.5,                   # pole powierzchni przekroju poprzecznego zbiornika
        'B': 0.25,                  # współczynnik wypływu
        'Q_o': [],                  # natężenie odpływu
    }
    fuzzy_val = {
        'k_e': 0.25,                # error
        'k_ce': 0.5,                # errorChange
        'k_u': 0.2,                 # u
    }
    N = 10000
    x = [0, ] # wzór na  sampletime to 1000 * sample time

    def JsonDate(self):
        listatoUaR = []
        with open('data.json') as json_data:
            data_dict = json.load(json_data)
            for key, value in data_dict.items():
                listatoUaR.append(value)
        listatoUaR.remove(listatoUaR[0])
        listatoUaR.remove(listatoUaR[-1])
        self.pid['sample_time'] = float(listatoUaR[0])
        self.pid['differential_time'] = float(listatoUaR[1]) # niewykorzystywane w fuzzy
        self.pid['integration_time'] = float(listatoUaR[2])  # niewykorzystywane w fuzzy
        self.pid['gain'] = float(listatoUaR[3])              # niewykorzystywane w fuzzy
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
        self.valvee['y'] = [[], [], [], [], []]
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
        self.JsonDate()
        for n in range(0, self.N):
            self.x.append(int(obiekt.pid['sample_time'] * n))
            self.pid_controler()
            self.valve()
            self.tankk()

    def display(self):
        self.count()
        plt.figure(1)
        plt.xlabel('t [s]')
        plt.ylabel('h [m]')
        plt.title('Wykres zależności poziomu substancji w zbiorniku - PID', fontsize=14)
        plt.plot(self.x, self.pid['h'], label='$PID$', color='blue')
        plt.axhline(self.pid['h_z'][-1], label='$Wartość _ zadana$', color='red', linewidth='2', linestyle='--')
        plt.legend()
        plt.savefig('static/plot.png')
        self.reset_data()


#FUZZY
    def Fuzzy(self):
        #fuzzyValue
        self.JsonDate()

        ins = ['DU', 'SU', 'MU', 'Z', 'MD', 'SD', 'DD']
        outs = ['BDU', 'DU', 'SU', 'MU', 'Z', 'MD', 'SD', 'DD', 'BDD']

        e = ctrl.Antecedent(np.arange(-1, 1, 0.01), 'e')
        ce = ctrl.Antecedent(np.arange(-1, 1, 0.01), 'ce')
        cu = ctrl.Consequent(np.arange(-1, 1, 0.01), 'cu')

        e['DU'] = fuzz.trimf(e.universe, [-1.333, -1, -0.6667])
        e['SU'] = fuzz.trimf(e.universe, [-1, -0.6667, -0.3333])
        e['MU'] = fuzz.trimf(e.universe, [-0.6667, -0.3333, -5.551e-17])
        e['Z'] = fuzz.trimf(e.universe, [-0.3333, 0.0, 0.3333])
        e['MD'] = fuzz.trimf(e.universe, [-5.551e-17, 0.3333, 0.6667])
        e['SD'] = fuzz.trimf(e.universe, [0.3333, 0.6667, 1])
        e['DD'] = fuzz.trimf(e.universe, [0.6667, 1, 1.333])

        ce['DU'] = fuzz.trimf(ce.universe, [-1.333, -1, -0.6667])
        ce['SU'] = fuzz.trimf(ce.universe, [-1, -0.6667, -0.3333])
        ce['MU'] = fuzz.trimf(ce.universe, [-0.6667, -0.3333, -5.551e-17])
        ce['Z'] = fuzz.trimf(ce.universe, [-0.3333, 0.0, 0.3333])
        ce['MD'] = fuzz.trimf(ce.universe, [-5.551e-17, 0.3333, 0.6667])
        ce['SD'] = fuzz.trimf(ce.universe, [0.3333, 0.6667, 1])
        ce['DD'] = fuzz.trimf(ce.universe, [0.6667, 1, 1.333])

        cu['BDU'] = fuzz.trimf(cu.universe, [-1.25, -1, -0.75])
        cu['DU'] = fuzz.trimf(cu.universe, [-1, -0.75, -0.5])
        cu['SU'] = fuzz.trimf(cu.universe, [-0.75, -0.5, -0.25])
        cu['MU'] = fuzz.trimf(cu.universe, [-0.5, -0.25, 0])
        cu['Z'] = fuzz.trimf(cu.universe, [-0.25, 0.0, 0.25])
        cu['MD'] = fuzz.trimf(cu.universe, [0, 0.25, 0.5])
        cu['SD'] = fuzz.trimf(cu.universe, [0.25, 0.5, 0.75])
        cu['DD'] = fuzz.trimf(cu.universe, [0.5, 0.75, 1])
        cu['BDD'] = fuzz.trimf(cu.universe, [0.75, 1, 1.25])

        # BDU
        rules = [ctrl.Rule(e[ins[0]] & ce[ins[0]], cu[outs[0]])]
        rules.append(ctrl.Rule(e[ins[1]] & ce[ins[0]], cu[outs[0]]))
        rules.append(ctrl.Rule(e[ins[0]] & ce[ins[1]], cu[outs[0]]))
        rules.append(ctrl.Rule(e[ins[1]] & ce[ins[1]], cu[outs[0]]))
        rules.append(ctrl.Rule(e[ins[2]] & ce[ins[0]], cu[outs[0]]))
        rules.append(ctrl.Rule(e[ins[0]] & ce[ins[2]], cu[outs[0]]))
        # DU
        rules.append(ctrl.Rule(e[ins[3]] & ce[ins[0]], cu[outs[1]]))
        rules.append(ctrl.Rule(e[ins[2]] & ce[ins[1]], cu[outs[1]]))
        rules.append(ctrl.Rule(e[ins[1]] & ce[ins[2]], cu[outs[1]]))
        rules.append(ctrl.Rule(e[ins[0]] & ce[ins[3]], cu[outs[1]]))
        # SU
        rules.append(ctrl.Rule(e[ins[0]] & ce[ins[4]], cu[outs[2]]))
        rules.append(ctrl.Rule(e[ins[1]] & ce[ins[3]], cu[outs[2]]))
        rules.append(ctrl.Rule(e[ins[2]] & ce[ins[2]], cu[outs[2]]))
        rules.append(ctrl.Rule(e[ins[3]] & ce[ins[1]], cu[outs[2]]))
        rules.append(ctrl.Rule(e[ins[4]] & ce[ins[0]], cu[outs[2]]))
        # MU
        rules.append(ctrl.Rule(e[ins[5]] & ce[ins[0]], cu[outs[3]]))
        rules.append(ctrl.Rule(e[ins[4]] & ce[ins[1]], cu[outs[3]]))
        rules.append(ctrl.Rule(e[ins[3]] & ce[ins[2]], cu[outs[3]]))
        rules.append(ctrl.Rule(e[ins[2]] & ce[ins[3]], cu[outs[3]]))
        rules.append(ctrl.Rule(e[ins[1]] & ce[ins[4]], cu[outs[3]]))
        rules.append(ctrl.Rule(e[ins[0]] & ce[ins[5]], cu[outs[3]]))
        # Z
        rules.append(ctrl.Rule(e[ins[0]] & ce[ins[6]], cu[outs[4]]))
        rules.append(ctrl.Rule(e[ins[1]] & ce[ins[5]], cu[outs[4]]))
        rules.append(ctrl.Rule(e[ins[2]] & ce[ins[4]], cu[outs[4]]))
        rules.append(ctrl.Rule(e[ins[3]] & ce[ins[3]], cu[outs[4]]))
        rules.append(ctrl.Rule(e[ins[4]] & ce[ins[2]], cu[outs[4]]))
        rules.append(ctrl.Rule(e[ins[5]] & ce[ins[1]], cu[outs[4]]))
        rules.append(ctrl.Rule(e[ins[6]] & ce[ins[0]], cu[outs[4]]))
        # MD
        rules.append(ctrl.Rule(e[ins[1]] & ce[ins[6]], cu[outs[5]]))
        rules.append(ctrl.Rule(e[ins[2]] & ce[ins[5]], cu[outs[5]]))
        rules.append(ctrl.Rule(e[ins[3]] & ce[ins[4]], cu[outs[5]]))
        rules.append(ctrl.Rule(e[ins[4]] & ce[ins[3]], cu[outs[5]]))
        rules.append(ctrl.Rule(e[ins[5]] & ce[ins[2]], cu[outs[5]]))
        rules.append(ctrl.Rule(e[ins[6]] & ce[ins[1]], cu[outs[5]]))
        # SD
        rules.append(ctrl.Rule(e[ins[6]] & ce[ins[2]], cu[outs[6]]))
        rules.append(ctrl.Rule(e[ins[5]] & ce[ins[3]], cu[outs[6]]))
        rules.append(ctrl.Rule(e[ins[4]] & ce[ins[4]], cu[outs[6]]))
        rules.append(ctrl.Rule(e[ins[3]] & ce[ins[5]], cu[outs[6]]))
        rules.append(ctrl.Rule(e[ins[2]] & ce[ins[6]], cu[outs[6]]))
        # DD
        rules.append(ctrl.Rule(e[ins[3]] & ce[ins[6]], cu[outs[7]]))
        rules.append(ctrl.Rule(e[ins[4]] & ce[ins[5]], cu[outs[7]]))
        rules.append(ctrl.Rule(e[ins[5]] & ce[ins[4]], cu[outs[7]]))
        rules.append(ctrl.Rule(e[ins[6]] & ce[ins[3]], cu[outs[7]]))
        # BDD
        rules.append(ctrl.Rule(e[ins[6]] & ce[ins[4]], cu[outs[8]]))
        rules.append(ctrl.Rule(e[ins[5]] & ce[ins[5]], cu[outs[8]]))
        rules.append(ctrl.Rule(e[ins[4]] & ce[ins[6]], cu[outs[8]]))
        rules.append(ctrl.Rule(e[ins[6]] & ce[ins[5]], cu[outs[8]]))
        rules.append(ctrl.Rule(e[ins[5]] & ce[ins[6]], cu[outs[8]]))
        rules.append(ctrl.Rule(e[ins[6]] & ce[ins[6]], cu[outs[8]]))

        cu_ctrl = ctrl.ControlSystem(rules)
        reg_cu = ctrl.ControlSystemSimulation(cu_ctrl)

        for x in range(0, self.N):
            self.pid['h_z'].append(obiekt.pid['h_z'][-1])

        for x in range(0, self.N):
            self.valvee['Q_d'].append(1)

        # h, Q_d, Q_o, h_z
        temp = [0]

        self.valvee['y'][0].append(0)
        self.valvee['y'][1].append(self.valvee['h_start'])
        self.valvee['y'][2].append(self.valvee['Q_d'][0])
        self.valvee['y'][3].append(-(self.tank['B']) * (sqrt(self.valvee['h_start'])))
        self.valvee['y'][4].append(self.pid['h_z'][0])
        e_delay = 0
        u_old = 0

        for x in range(1, self.N):
            e = self.pid['h_z'][x] - self.valvee['y'][1][x - 1]
            ke = e * (1 / self.fuzzy_val['k_e'])
            temp.append(ke)
            if e > 1:
                e = 1
            elif e < -1:
                e = -1

            ce = e - e_delay
            ce *= (1 / obiekt.pid['sample_time'])  # T_p lub T_0
            ce *= (1 / self.fuzzy_val['k_ce'])
            if ce > 1:
                ce = 1
            elif ce < -1:
                ce = -1

            reg_cu.input['e'] = e # error
            reg_cu.input['ce'] = ce # errorChange
            reg_cu.compute()
            cu = reg_cu.output['cu']
            u = (cu * obiekt.pid['sample_time'] + u_old)  # T_p lub T_0
            temp.append(u)

            self.valvee['Q_d'][x] = self.fuzzy_val['k_u'] * u

            if (self.valvee['Q_d'][x] >self.valvee['Q_d_max']):
                self.valvee['Q_d'][x] = self.valvee['Q_d_max']
            elif (self.valvee['Q_d'][x] < self.valvee['Q_d_min']):
                self.valvee['Q_d'][x] = self.valvee['Q_d_min']

            self.valvee['y'][0].append((x) * obiekt.pid['sample_time'])
            self.valvee['y'][2].append(self.valvee['Q_d'][x])
            self.valvee['y'][3].append(self.tank['B'] * (sqrt(self.valvee['y'][1][x - 1])))
            self.valvee['y'][1].append(((self.valvee['y'][2][x] - self.valvee['y'][3][x]) * obiekt.pid['sample_time']) / self.tank['A'] + self.valvee['y'][1][x - 1])
            self.valvee['y'][4].append(self.pid['h_z'][x])

            if self.valvee['y'][1][len(self.valvee['y'][1]) - 1] > self.tank['h_max']:
                return('Przekroczono górny limit', 'Ostrzeżenie', 'warning')

            if self.valvee['y'][1][len(self.valvee['y'][1]) - 1] < self.tank['h_min']:
                return('Przekroczono dolny limit', 'Ostrzeżenie', 'warning')
            e_delay = e
            u_old = u

    def FuzzyDisplay(self):
        self.Fuzzy()
        plt.figure(2)
        plt.xlabel('t [s]')
        plt.ylabel('h [m]')
        plt.title('Wykres zależności poziomu substancji w zbiorniku - PID Fuzzy', fontsize=14)
        plt.plot(self.valvee['y'][0], self.valvee['y'][1], label='$fuzzyPID$', color='green')
        plt.axhline(self.pid['h_z'][-1], label='$Wartość _ zadana$', color='red', linewidth='2', linestyle='--')
        plt.legend()
        plt.savefig('static/Fuzze.png')
        self.reset_data()

    def Porownanie(self):
        self.count()
        self.Fuzzy()
        plt.figure(2)
        plt.xlabel('t [s]')
        plt.ylabel('h [m]')
        plt.title('PID vs Fuzzy PID', fontsize=14)
        plt.plot(self.x, self.pid['h'], label='$PID$', color='blue')
        plt.plot(self.valvee['y'][0], self.valvee['y'][1], label='$fuzzyPID$', color='green')
        plt.axhline(self.pid['h_z'][-1], label='$Wartość _ zadana$', color='red', linewidth='2', linestyle='--')
        plt.legend()
        plt.savefig('static/porownanie.png')
        self.reset_data()

obiekt = UAR()
UAR.Porownanie(obiekt)