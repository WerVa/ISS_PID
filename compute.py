from math import sqrt
import matplotlib.pyplot as plt
import json


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
        self.pid['error'] = [0, ]
        self.pid['gain'] = 1.0
        self.pid['sample_time'] = 0.05
        self.pid['differential_time'] = 0.05
        self.pid['integration_time'] = 0.75
        self.pid['h_z'] = []
        self.pid['h'] = [0, ]

        self.valvee['u_max'] = 10
        self.valvee['u_min'] = -10
        self.valvee['u'] = [0.1, ]
        self.valvee['Q_d_min'] = 0
        self.valvee['Q_d_max'] = 1
        self.valvee['Q_d'] = []

        self.tank['h_max'] = 10
        self.tank['h_min'] = -10
        self.tank['A'] = 2.5
        self.tank['B'] = 0.25
        self.tank['Q_o'] = []

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


obiekt = UAR()