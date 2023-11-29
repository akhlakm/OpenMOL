""" OpenMOL utility functions and classes.

	This file is a part of OpenMOL python module.
	License GPLv3.0 Copyright (c) 2023 Akhlak Mahmood """

import re

class AttrDict(dict):
    """ Adds a convenient way to access dictionary items
    as properties """

    def __init__(self, d : dict = {}):
        if type(d) == dict:
            for k,v in d.items():
                self.__setattr__(k, v)
        else:
            raise ValueError("Initial value must be a dict")

    def __getattr__(self, key : str) -> any:
        return self[key]

    def __setattr__(self, key : str, value):
        self.__setitem__(key, value)

    def __setitem__(self, key : str, value):
        # do not allow any key from dict's namespace
        if key in dir({}):
            raise KeyError(key)

        # key has to be a string
        if type(key) != str:
            raise KeyError(key)

        # key can't start with a digit, must be alphanumeric
        # unscore allowed
        search = re.compile(r'^[a-zA-Z_][a-zA-Z_0-9]*$').search
        if not bool(search(key)):
            raise KeyError(key)

        super(AttrDict, self).__setitem__(key, value)

    def __dir__(self):
        return super().__dir__() + [str(k) for k in self.keys()]

