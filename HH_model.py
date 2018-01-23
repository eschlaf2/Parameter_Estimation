import models
import numpy as np
from integrators import ruku4

default_params = {
    'x0': -1.6,  # [s]
    'y0': 1.,
    'tau0': 2857,
    'tau1': 1.0,
    'tau2': 10.,
    'Irest1': 3.1,
    'Irest2': 0.45,
    'gamma': 1e-2,
    'x1_init': 0.,
    'y1_init': -5.,
    'z_init': 3.,
    'x2_init': 0.,
    'y2_init': 0.,
    'g_init': 0.,
    'observation_sigmas': None,
    'noise_ensemble1': 25e-3,
    'noise_ensemble2': 25e-2,
    'a': 5.,
    'b': 4.,
    'c': 0.3,
    'd': 3.5}


def f1(self, x1, x2, z):
    k = 3
    return (x1**3 - k * x1**2) * (x1 < 0) + \
        (x1 * (x2 - 0.6 * (z - 4)**2)) * (x1 >= 0)  # k=3


def f2(self, x2):
    return 0. * (x2 < -0.25) + (6 * (x2 + 0.25)) * (x2 >= -0.25)


class HH_model(models.Model):
    def __init__(self, params=None,
                 total_time=2500, dt_sample=0.1, dt_integrate=None,
                 **kwargs):
        if not params:
            params = {}
            for key, value in default_params.iteritems():
                params[key] = value
        for key, value in kwargs.iteritems():
            params[key] = float(value)
        self.params = params

        self.initial_conditions = [params['x1_init'], params['y1_init'],
                                   params['z_init'], params['x2_init'],
                                   params['y2_init'], params['g_init']]
        self.noise = [params['noise_ensemble1'],
                      0.0,
                      0.,
                      params['noise_ensemble2'],
                      0.,
                      0.]

        self.integrator = 'ruku4'

        self.var_names = ['x1', 'y1', 'z', 'x2', 'y2', 'g']
        self.parameter_names = ['$I_{ext1}$', '$I_{ext2}$', '$I_{extz}$']

        self.num_params = len(self.parameter_names)
        self.num_state_vars = len(self.var_names)
        self.num_observations = 1
        self.num_augmented_states = self.num_params + self.num_state_vars

        # self._total_time = total_time  # 2500 for epileptor, 160 for FN
        # self._dt_sample = dt_sample
        # self._dt_integrate = self._dt_sample if dt_integrate is None \
        #     else dt_integrate
        # self._num_samples = int(models.set_num_samples(
        #     self._total_time, self._dt_sample))
        # self._steps_per_sample = int(models.set_steps_per_sample(
        #     self._dt_sample, self._dt_integrate))

        I_ext1 = (
            self.Irest1 * np.ones(self._num_samples)).reshape(
            1, self._num_samples)
        I_ext2 = (
            self.Irest2 * np.ones(self._num_samples)).reshape(
            1, self._num_samples)
        I_extz = (
            1. * np.ones(self._num_samples)).reshape(
            1, self._num_samples)
        self.parameters = np.vstack((I_ext1, I_ext2, I_extz))

    def dynamics(self, state, time_varying_params):
        '''
        TODO: unscented kalman filter parameter
        '''
        x1, y1, z, x2, y2, g = state
        I_ext1, I_ext2, I_extz = time_varying_params

        x1_dot = y1 - f1(x1, x2, z) - z + I_ext1  # self.Irest1
        y1_dot = (self.y0 - self.a * x1 * x1 - y1) / \
            self.tau1  # a = 5., tvb param d
        z_dot = I_extz / self.tau0 * \
            (self.b * (x1 - self.x0) - z)  # b = 4., tvb const
        x2_dot = -y2 + x2 - x2**3 + I_ext2 + \
            2. * g - self.c * (z - self.d)  # + self.Irest2 c = 0.3, d = 3.5
        y2_dot = (-y2 + f2(x2)) / self.tau2
        g_dot = -self.gamma * (g - 0.1 * x1)

        return np.array([x1_dot, y1_dot, z_dot, x2_dot, y2_dot, g_dot])

    def observation_function(self, augmented_state):
        x1 = augmented_state[self.num_params, :]
        x2 = augmented_state[self.num_params + 3, :]
        return -x1 + x2

    def set_initial_estimate(self, initial_estimate):
        x1 = self.noisy_data[0, 0] / 2.
        x2 = self.noisy_data[0, 0] / 2. + x1
        initial_estimate[self.num_params] = x1
        initial_estimate[self.num_params + 3] = x2
        return initial_estimate

    def transition_function(self, augmented_state):
        parameters, state = np.split(augmented_state, [self.num_params, ])
        state = ruku4(self.dynamics,
                      state,
                      parameters,
                      self.dt_integrate,
                      self.steps_per_sample,
                      self.noise)
        return np.vstack((parameters, state))
