# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 17:33:41 2013.

Utility functions for calcium-based STDP using simplified calcium model as in
(Graupner and Brunel, 2012).

@author: Giuseppe Chindemi
@remark: Copyright © BBP/EPFL 2005-2016; All rights reserved.
         Do not distribute without further notice.
"""

# pylint: disable=R0914, R0912

import logging
import numpy as np
import integrators
import pickle
import matplotlib.pyplot as plt
from integrators import ruku4
from pyedflib import EdfReader

logging.basicConfig(level=logging.WARN)

# Note: having debug logging statements increases the run time by ~ 25%,
# because they exist in tight loops, and expand their outputs, even when
# debug is off, so we disable logging if possible.  Set this to true if
# verbose output is needed
LOGGING_DEBUG = False

plot = True


def logging_debug_vec(fmt, vec):
    '''log to debug a vector'''
    if LOGGING_DEBUG:
        logging.debug(fmt, ', '.join(map(str, vec)))


def logging_debug(*args):
    '''wrapper to log to debug a vector'''
    if LOGGING_DEBUG:
        print(args)
    # if LOGGING_DEBUG:
    #     logging.debug(*args)


def load_obj(name):
    with open(name, 'rb') as f:
        return pickle.load(f)


def read_edf_file(filename):
    f = EdfReader(filename)
    chan = 0
    data = f.readSignal(chan) / 1000.
    sample_freq = f.getSampleFrequency(chan)
    f._close()
    return data, sample_freq


def read_pkl_file(filename):
    f = load_obj(filename)
    return f.noisy_data[0].reshape(-1), int(round(1. / f.dt_sample))


def subsample_data(data, sample_freq, dt_sample):
    dt_sample = max(dt_sample, 1. / sample_freq)
    dstep = int(sample_freq * dt_sample)
    data = data[::dstep]
    num_samples = len(data)
    return data, dt_sample, num_samples


def read_file(filename):
    filetype_switcher = {
        'edf': read_edf_file,
        'pkl': read_pkl_file
    }
    filetype = filename.split('.')[-1]
    data, sample_freq = filetype_switcher[filetype](filename)
    return data, sample_freq


class Protocol(object):

    """Protocol"""

    def __init__(self, params=None,
                 total_time=2500, prot_id=None, **kwargs):
        self.prot_id = prot_id
        self.params = dict(params)
        self.total_time = total_time
        for key, value in kwargs.iteritems():
            self.params[key] = float(value)


class Model(object):
    def __init__(self):
        self.augmented_state = []

    def _get_total_time(self):
        return self._total_time

    def _set_total_time(self, value):
        self._total_time = value
        self._num_samples = set_num_samples(
            self._total_time, self._dt_sample)

    def _get_dt_sample(self):
        return self._dt_sample

    def _set_dt_sample(self, value):
        self._dt_sample = value
        self._steps_per_sample = \
            set_steps_per_sample(self._dt_sample, self._dt_integrate)
        self._num_samples = \
            set_num_samples(self._total_time, self._dt_sample)

    def _get_dt_integrate(self):
        return self._dt_integrate

    def _set_dt_integrate(self, value):
        self._dt_integrate = value
        self._steps_per_sample = \
            set_steps_per_sample(self._dt_sample, self._dt_integrate)

    total_time = property(_get_total_time, _set_total_time)
    dt_sample = property(_get_dt_sample, _set_dt_sample)
    dt_integrate = property(_get_dt_integrate, _set_dt_integrate)
    num_samples = property(lambda self: self._num_samples)
    steps_per_sample = property(lambda self: self._steps_per_sample)
    time = property(lambda self:
                    np.arange(0, self._total_time, self._dt_sample))

    def generate_simulation(self, plot=plot):
        '''
            Simulates true and noisy trajectory based on previously
            defined model and parameter functions
            (Uses global vars)
            '''

        # Simulate model
        true_state = self.integrate_model()
        self.augmented_state = np.vstack((self.parameters * np.ones(
            (self.dims_params, self._num_samples)), true_state))\
            if self.dims_params > 0 \
            else true_state
        self.dims_augmented_state = self.dims_params + self.dims_state_vars

        # Observation noise
        if self.observation_sigmas is None:
            self.observation_sigmas = list(np.sqrt([0.2 * 0.2 * np.var(
                self.observation_function(self.augmented_state))]))
        if type(self.observation_sigmas) != list:
            self.observation_sigmas = [self.observation_sigmas]
        if len(self.observation_sigmas) > self.dims_observations:
            self.observation_sigmas = \
                self.observation_sigmas[:self.dims_observations]
        while len(self.observation_sigmas) < self.dims_observations:
            self.observation_sigmas.append(0.)

        observation_noise = np.diag(self.observation_sigmas)

        # Create noisy data from true trajectory
        self.data = self.observation_function(self.augmented_state)
        self.noisy_data = self.data + \
            np.matmul(observation_noise,
                      np.random.randn(self.dims_observations,
                                      self._num_samples))

        if plot:
            self.plot_simulated_data()
        return self.noisy_data  # if noisy else self.data

    def integrate(self, state, time_varying_params):
        switcher = {
            'ruku4': integrators.ruku4,
            'euler': integrators.euler,
            'euler_maruyama': integrators.euler_maruyama,
            'test_integrator': integrators.test_integrator
        }
        logging_debug('state.shape: ', state.shape)
        logging_debug('time_varying_params.shape: ', time_varying_params.shape)
        return switcher[self.integrator](self.model_function,
                                         state,
                                         time_varying_params,
                                         self.dt_integrate,
                                         self.steps_per_sample,
                                         self.noise)

    def integrate_model(self):
        true_state = np.zeros((len(self.initial_conditions),
                               self._num_samples))  # allocate
        true_state[:, 0] = self.initial_conditions
        for n in range(self._num_samples - 1):
            logging_debug('n: ', n)
            x_temp = true_state[:, n]
            true_state[:, n + 1] = \
                self.integrate(state=x_temp,
                               time_varying_params=self.parameters[:, n])
            if any(np.isinf(true_state[:, n + 1])):
                break
        return true_state

    def plot_simulated_data(self):
        '''Plot simulation'''
        plt.rc('text', usetex=True)
        plt.figure(figsize=(10, 2))
        plt.plot(self.time, self.noisy_data[0, :],
                 # plt.plot(self.target[0],
                 'bd', markeredgecolor='blue',
                 mfc='blue', ms=3, label='noisy data')
        plt.plot(self.time, self.observation_function(self.augmented_state).T,
                 'k', linewidth=2, label='actual')
        # plt.figure()
        # for i in range(self.dims_state_vars):
        # plt.plot(self.augmented_state[self.dims_params + 2, :],
        #          label=self.var_names[2])
        plt.xlabel('t')
        plt.legend()
        plt.axis('tight')
        plt.title('Simulation')


def set_num_samples(total_time, dt_sample):
    return round(total_time / dt_sample)


def set_steps_per_sample(dt_sample, dt_integrate):
    return round(dt_sample / dt_integrate)


def load_protocols(filename=None, plot=plot, total_time=None, dt_sample=0.1):
    """Load simulation using params from Jirsa, 2014."""

    # protocols = [Protocol(prot_id='default'),
    #              Protocol(prot_id='clean', observation_sigmas=0., tau0=1000),
    #              Protocol(prot_id='noiseless',
    #                       observation_sigmas=0.,
    #                       noise_ensemble1=0.,
    #                       noise_ensemble2=0.,
    #                       tau0=2000)]
    # target = [epileptor_model(params=protocols[0].params,
    #                           total_time=total_time, dt_sample=dt_sample).
    #           generate_simulation(plot=plot),
    #           epileptor_model(params=protocols[1].params,
    #                           total_time=total_time, dt_sample=dt_sample).
    #           generate_simulation(plot=plot),
    #           epileptor_model(params=protocols[2].params,
    #                           total_time=total_time, dt_sample=dt_sample).
    #           generate_simulation(plot=plot)]

    prot_id = filename.split('/')[-1]
    data, sample_freq = read_file(filename)
    data, dt_sample, num_samples = \
        subsample_data(data, sample_freq, dt_sample)
    if total_time is None:
        total_time = int(num_samples * dt_sample)
        print(total_time)
    else:
        total_time = int(min(total_time, num_samples // dt_sample))

    protocols = [Protocol(prot_id=prot_id, total_time=total_time)]
    target = [data[:int(total_time / dt_sample)]]

    # f = EdfReader(
    #     '/Users/emilyschlafly/BU/Kramer_rotation/ieeg_data/' +
    #     'I002_A0003_D010/outputEdf_EDF/outputEdf_0.edf')
    # chan = 0
    # data = f.readSignal(chan)
    # sample_freq = f.getSampleFrequency(chan)
    # if not total_time:
    #     num_samples = len(data)
    #     total_time = num_samples / sample_freq
    # data = data[::int(round(sample_freq * dt_sample))]
    # data = data[:int(total_time / dt_sample)] / 1000.
    # total_time = min(total_time, len(data) * dt_sample)

    # f._close()
    # del f

    # if plot:
    #     plt.plot(np.arange(0, total_time, dt_sample), data)

    # target = [data]
    # protocols = [Protocol(prot_id='1', total_time=total_time)]

    return protocols, target


def rmse(estimate, target):
    logging_debug('estimate: %f, %s', estimate, type(estimate))
    logging_debug('target: %f, %s', target, type(target))
    return np.sqrt(((estimate - target)**2).mean())


def protocol_outcome(protocol, param=param_epileptor, plot=plot):
    """Compute the average synaptic gain for a given stimulation protocol and
    model parameters.

    :param protocol: epileptor_util.Protocol
        The stimulation protocol.
    :param model: dict
        Parameters of the Epileptor model
    """

    estimate = epileptor_model(param, total_time=protocol.total_time).\
        generate_simulation(plot=plot)

    return estimate

# ================================================
# EXAMPLE OF MODEL (Fitzhugh-Nagumo)
# ================================================
class FN_model(Model):
    def __init__(self, params=None,
                 total_time=160, dt_sample=0.2, dt_integrate=None,
                 **kwargs):
        if not params:
            params = {
                'a': 0.7,
                'b': 0.8,
                'c': 3.,
                'observation_sigmas': 25e-2
            }
        for key, value in kwargs.iteritems():
            params[key] = float(value)
        self.params = params
        self.a, self.b, self.c = params['a'], params['b'], params['c']

        self.initial_conditions = [0., 0.]
        self.noise = [0., 0.]
        self.observation_sigmas = params['observation_sigmas']

        self.integrator = 'ruku4'

        self.var_names = ['v', 'w']
        # self.parameter_names = ['$I_{ext}$', 'a', 'b', 'c']
        self.parameter_names = ['$I_{ext}$', 'a']
        # self.parameter_names = ['$I_{ext}$']

        self.dims_params = len(self.parameter_names)
        self.dims_state_vars = len(self.var_names)
        self.dims_observations = 1
        self.dims_augmented_state = self.dims_params + self.dims_state_vars

        self._total_time = total_time  # 2500 for epileptor, 160 for FN
        self._dt_sample = dt_sample
        self._dt_integrate = self._dt_sample if dt_integrate is None \
            else dt_integrate
        self._num_samples = int(set_num_samples(
            self._total_time, self._dt_sample))
        self._steps_per_sample = int(set_steps_per_sample(
            self._dt_sample, self._dt_integrate))

        Iext = np.arange(1, self._num_samples + 1) / 250. * 2 * np.pi
        Iext = -0.4 - 1.01 * (np.abs(np.sin(Iext / 2.)))
        a = (
            self.a * np.ones(self._num_samples)).reshape(
            1, self._num_samples)
        self.parameters = np.vstack((Iext, a))

    def set_initial_estimate(self, initial_estimate):
        w = self.noisy_data[0, 0]
        initial_estimate[-2] = w
        return initial_estimate

    def model_function(self, state, parameters):
        a = self.a
        b = self.b
        c = self.c
        v, w = state
        # input_current = parameters[0]
        input_current, a = parameters
        logging_debug('model_function -> state: ', state)
        logging_debug('model_function -> v: ', v)
        logging_debug('model_function -> v.shape', v.shape)
        logging_debug('model_function -> input_current.shape',
                      input_current.shape)
        # input_current, a, b, c = parameters
        v_dot = c * (w + v - v**3 / 3 + input_current)
        w_dot = -(v - a + b * w) / c
        return np.array([v_dot, w_dot])

    def observation_function(self, augmented_state):
        w = augmented_state[-2, :]
        return w

    def transition_function(self, augmented_state):
        parameters, state = np.split(augmented_state, [self.dims_params, ])
        state = ruku4(self.model_function,
                      state,
                      parameters,
                      self.dt_integrate,
                      self.steps_per_sample,
                      self.noise)
        return np.vstack((parameters, state))

