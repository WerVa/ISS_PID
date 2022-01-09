import matplotlib.pyplot as plt
import math
import numpy as np
import skfuzzy as fuzz
from skfuzzy import control as ctrl
h_z_const = 5;  #OK

beta = 0.25;    #OK
A = 2.5;        #OK
h_0 = 0;        #OK
h_min= 0;       #OK
h_max = 10;     #OK

Q_d_min = 0;    #OK
Q_d_max = 1;    #OK

T_p = 0.05;     #OK
simtime = 500;  #OK

k_e = 0.25;
T_0 = 0.01;
k_ce = 0.5;
k_u = 0.1;      #OK

ins = ['DU', 'SU', 'MU', 'Z', 'MD', 'SD', 'DD']
outs = ['BDU','DU', 'SU', 'MU', 'Z', 'MD', 'SD', 'DD', 'BDD']

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
#DU
rules.append(ctrl.Rule(e[ins[3]] & ce[ins[0]], cu[outs[1]]))
rules.append(ctrl.Rule(e[ins[2]] & ce[ins[1]], cu[outs[1]]))
rules.append(ctrl.Rule(e[ins[1]] & ce[ins[2]], cu[outs[1]]))
rules.append(ctrl.Rule(e[ins[0]] & ce[ins[3]], cu[outs[1]]))
#SU
rules.append(ctrl.Rule(e[ins[0]] & ce[ins[4]], cu[outs[2]]))
rules.append(ctrl.Rule(e[ins[1]] & ce[ins[3]], cu[outs[2]]))
rules.append(ctrl.Rule(e[ins[2]] & ce[ins[2]], cu[outs[2]]))
rules.append(ctrl.Rule(e[ins[3]] & ce[ins[1]], cu[outs[2]]))
rules.append(ctrl.Rule(e[ins[4]] & ce[ins[0]], cu[outs[2]]))
#MU
rules.append(ctrl.Rule(e[ins[5]] & ce[ins[0]], cu[outs[3]]))
rules.append(ctrl.Rule(e[ins[4]] & ce[ins[1]], cu[outs[3]]))
rules.append(ctrl.Rule(e[ins[3]] & ce[ins[2]], cu[outs[3]]))
rules.append(ctrl.Rule(e[ins[2]] & ce[ins[3]], cu[outs[3]]))
rules.append(ctrl.Rule(e[ins[1]] & ce[ins[4]], cu[outs[3]]))
rules.append(ctrl.Rule(e[ins[0]] & ce[ins[5]], cu[outs[3]]))
#Z
rules.append(ctrl.Rule(e[ins[0]] & ce[ins[6]], cu[outs[4]]))
rules.append(ctrl.Rule(e[ins[1]] & ce[ins[5]], cu[outs[4]]))
rules.append(ctrl.Rule(e[ins[2]] & ce[ins[4]], cu[outs[4]]))
rules.append(ctrl.Rule(e[ins[3]] & ce[ins[3]], cu[outs[4]]))
rules.append(ctrl.Rule(e[ins[4]] & ce[ins[2]], cu[outs[4]]))
rules.append(ctrl.Rule(e[ins[5]] & ce[ins[1]], cu[outs[4]]))
rules.append(ctrl.Rule(e[ins[6]] & ce[ins[0]], cu[outs[4]]))
#MD
rules.append(ctrl.Rule(e[ins[1]] & ce[ins[6]], cu[outs[5]]))
rules.append(ctrl.Rule(e[ins[2]] & ce[ins[5]], cu[outs[5]]))
rules.append(ctrl.Rule(e[ins[3]] & ce[ins[4]], cu[outs[5]]))
rules.append(ctrl.Rule(e[ins[4]] & ce[ins[3]], cu[outs[5]]))
rules.append(ctrl.Rule(e[ins[5]] & ce[ins[2]], cu[outs[5]]))
rules.append(ctrl.Rule(e[ins[6]] & ce[ins[1]], cu[outs[5]]))
#SD
rules.append(ctrl.Rule(e[ins[6]] & ce[ins[2]], cu[outs[6]]))
rules.append(ctrl.Rule(e[ins[5]] & ce[ins[3]], cu[outs[6]]))
rules.append(ctrl.Rule(e[ins[4]] & ce[ins[4]], cu[outs[6]]))
rules.append(ctrl.Rule(e[ins[3]] & ce[ins[5]], cu[outs[6]]))
rules.append(ctrl.Rule(e[ins[2]] & ce[ins[6]], cu[outs[6]]))
#DD
rules.append(ctrl.Rule(e[ins[3]] & ce[ins[6]], cu[outs[7]]))
rules.append(ctrl.Rule(e[ins[4]] & ce[ins[5]], cu[outs[7]]))
rules.append(ctrl.Rule(e[ins[5]] & ce[ins[4]], cu[outs[7]]))
rules.append(ctrl.Rule(e[ins[6]] & ce[ins[3]], cu[outs[7]]))
#BDD
rules.append(ctrl.Rule(e[ins[6]] & ce[ins[4]], cu[outs[8]]))
rules.append(ctrl.Rule(e[ins[5]] & ce[ins[5]], cu[outs[8]]))
rules.append(ctrl.Rule(e[ins[4]] & ce[ins[6]], cu[outs[8]]))
rules.append(ctrl.Rule(e[ins[6]] & ce[ins[5]], cu[outs[8]]))
rules.append(ctrl.Rule(e[ins[5]] & ce[ins[6]], cu[outs[8]]))
rules.append(ctrl.Rule(e[ins[6]] & ce[ins[6]], cu[outs[8]]))

cu_ctrl = ctrl.ControlSystem(rules)
reg_cu = ctrl.ControlSystemSimulation(cu_ctrl)

samples = (int)(simtime/T_p)+1;
h_z = [];
for x in range(samples):
    h_z.append(h_z_const);
    
Q_d = [];
for x in range(samples):
    Q_d.append(1);

y = [[],[],[],[], []];
#h, Q_d, Q_o, h_z
temp = [0]

y[0].append(0);
y[1].append(h_0);
y[2].append(Q_d[0]);
y[3].append(beta*(math.sqrt(h_0)));
y[4].append(h_z[0])
e_delay = 0
u_old = 0


for x in range(1,samples): 

    e = h_z[x] - y[1][x-1]
    ke = e*(1/k_e)
    temp.append(ke);
    if (e > 1):
        e = 1
    elif(e < -1):
        e = -1

    ce = e-e_delay    #MAM
    ce *= (1/T_0) #T_p lub T_0
    ce *= (1/k_ce)
    if (ce > 1):
        ce = 1
    elif(ce < -1):
        ce = -1

    
    reg_cu.input['e'] = e
    reg_cu.input['ce'] = ce

    reg_cu.compute()

    cu = reg_cu.output['cu']
    
    u = (cu*T_0 + u_old)  #T_p lub T_0

    #temp.append(u);
    
    Q_d[x] = k_u*u
    
    if (Q_d[x] > Q_d_max):
        Q_d[x] = Q_d_max
    elif(Q_d[x] < Q_d_min):
        Q_d[x] = Q_d_min
    
    y[0].append((x)*T_p);
    y[2].append(Q_d[x]);
    y[3].append(beta*(math.sqrt(y[1][x-1])));
    y[1].append(((y[2][x]-y[3][x])*T_p)/A + y[1][x-1]);
    y[4].append(h_z[x])
    if  y[1][len(y[1])-1] > h_max:
       print('Przekroczono górny limit','Ostrzeżenie','warning');
       break;
    
    if  y[1][len(y[1])-1] < h_min:
       print('Przekroczono dolny limit','Ostrzeżenie','warning')
       break;
    e_delay = e
    u_old = u

fig = plt.figure;

plt.title('Wykres zależności poziomu substancji w zbiorniku', fontsize=14)
plt.plot(y[0], y[1], label='$h$', color='blue')
plt.plot(y[0], y[4], label='$h_z$', color='red', linewidth='2', linestyle='--')



plt.show();


