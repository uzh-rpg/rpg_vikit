#!/usr/bin/python

import os

class RosNode:
  def __init__(self, package, executable):
    self._package = package
    self._executable = executable
    
  def run(self, parameter_dictionary):
    parameter_string = ''
    for key in parameter_dictionary.keys():
      parameter_string = parameter_string + ' _' + key + ':=' + str(parameter_dictionary[key])
    print('Starting ROS node.')
    os.system('rosrun ' + self._package + ' ' + self._executable + ' ' + parameter_string)
    print('ROS node finished processing.')
