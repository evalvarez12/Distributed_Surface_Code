import qutip as qt
import tools.operations


def generate_noisy_pair(F):
    a = qt.bell_state('00') * qt.bell_state('00').dag()
    b = qt.bell_state('01') * qt.bell_state('01').dag() \
        + qt.bell_state('10') * qt.bell_state('10').dag() \
        + qt.bell_state('11') * qt.bell_state('11').dag()
    W = F*a + (1 - F)/3.*b
    H = qt.snot(2, 1)
    return H*W*H.dag()


def generate_werner(F):
    a = qt.bell_state('00')
    a = a * a.dag()
    return F*a + (1-F)/4.*qt.qeye(4)


Fi = .9
State_initial = generate_noisy_pair(Fi)
State = qt.tensor(State_initial, State_initial)
CNOT = qt.cnot(4, 0, 2) * qt.cnot(4, 3, 1)
State = CNOT * State * CNOT.dag()

H = qt.snot(4, 1) * qt.snot(4, 3)
State = H * State * H.dag()

m1, State_collapsed1 = operations.measure_single(State, 4, 2)
m2, State_collapsed = operations.measure_single(State_collapsed1, 3, 2)
print(m1, m2)

State_final = qt.snot(2, 1) * State_collapsed * qt.snot(2, 1).dag()


State_ref = qt.snot(2, 1) * qt.bell_state('00')
State_ref = State_ref * State_ref.dag()

# initial Fidelity
F_initial = (State_ref * State_initial).tr()
F_final = (State_ref * State_final).tr()


def single_selection(F1, F2):
    F = F1*F2 + (1 - F1)*(1 - F2)/9
    F = F/(F1*F2 + F1*(1 - F2)/3 + F2*(1 - F1)/3 + 5*(1 - F1)*(1 - F2)/9)
    return F


print(single_selection(Fi, Fi), F_final)
