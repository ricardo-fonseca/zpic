#!/usr/bin/env python3

import math
import numpy as np
import scipy.special

nm = 1.e-9
um = 1.e-6
mm = 1.e-3
cm = 1.e-2
m = 1.
uA = 1.e-6
mA = 1.e-3
A = 1.
km = 1.e3
nT = 1.e-9
uT = 1.e-6
mT = 1.e-3
T = 1.

def normalize(v):
  l = np.linalg.norm(v)
  epsilon = 1e-9
  if l < epsilon:
    raise ValueError("Vector with |v| < e normalized")
  return v / l

# filamentary current loop
class Loop(object):
  def __init__(self, position, normal, radius, current):
    self.p = np.array(position, dtype = np.float64)
    self.n = normalize(np.array(normal, dtype = np.float64))
    self.r = np.float64(radius)
    self.i = np.float64(current)

# vector magnetic field (B) produced by collection of current loops
class Field(object):
  def __init__(self,
               length_units = m,
               current_units = A,
               field_units = T):
    self.loops = []
    self._length_units = length_units
    self._field_units = field_units
    self._epsilon = np.finfo(np.float64).eps
    
  def addLoop(self, loop):
    self.loops.append(loop)
    
  def evaluate(self, position):
    _p = np.atleast_2d(position)
    B = np.zeros(_p.shape)
    for loop in self.loops:
      B += self._evalLoop(_p, loop)
    return np.squeeze(B / self._field_units)
  
  def _evalLoop(self, p, loop):
    r_vect = (p - loop.p) * self._length_units
    r = np.linalg.norm(r_vect, axis=1, keepdims=True)
    z = r_vect.dot(loop.n.T)
    rho_vect = r_vect - np.outer(z, loop.n)
    rho = np.linalg.norm(rho_vect, axis=1)
    rho_vect[rho > self._epsilon,] = \
      (rho_vect[rho > self._epsilon,].T/rho[rho > self._epsilon]).T
    
    a = loop.r * self._length_units
    alpha2 = a*a + rho*rho + z*z - 2.*a*rho
    beta2 = a*a + rho*rho + z*z + 2.*a*rho
    beta = np.sqrt(beta2)        
    c = 4.e-7 * loop.i  # \mu_0  I / \pi
    a2b2 = alpha2 / beta2
    Ek2 = scipy.special.ellipe(1. - a2b2)
    Kk2 = scipy.special.ellipkm1(a2b2)
    
    denom = (2. * alpha2 * beta * rho)
    with np.errstate(invalid='ignore'):
      numer = c*z*((a*a + rho*rho + z*z)*Ek2 - alpha2*Kk2)
    sw = np.abs(denom) > self._epsilon
    Brho = np.zeros(numer.shape)
    Brho[sw] = numer[sw] / denom[sw]

    denom = (2. * alpha2 * beta)
    with np.errstate(invalid='ignore'):
      numer = c*((a*a - rho*rho - z*z)*Ek2 + alpha2*Kk2)
    sw = np.abs(denom) > self._epsilon
    Bz = np.full(numer.shape, np.inf)
    Bz[sw] = numer[sw] / denom[sw]

    return (Brho * rho_vect.T).T + np.outer(Bz, loop.n)
