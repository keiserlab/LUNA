from pymol import cmd
from pymol import util

from util.file import (get_filename, get_file_format)

import logging

logger = logging.getLogger()


class PymolWrapper:

    def load(self, input_file, obj_name=None):
        self.input_file = input_file
        if not obj_name:
            obj_name = get_filename(input_file)
        cmd.load(input_file, obj_name)

    def show(self, tuples):
        for representation, selection in tuples:
            cmd.show(representation, selection)

    def hide(self, tuples):
        for representation, selection in tuples:
            cmd.hide(representation, selection)

    def hide_all(self):
        self.hide([('everything', '')])

    def add_pseudoatom(self, name, opts=None):
        opts = opts or {}
        cmd.pseudoatom(name, **opts)

    def distance(self, name, sel1, sel2):
        cmd.distance(name, sel1, sel2)

    def save_png(self, output_file, width=1200, height=1200, dpi=100, ray=1):
        cmd.png(output_file, width, height, dpi, ray)

    def save_session(self, output_file):
        if get_file_format(output_file) != 'pse':
            output_file += '.pse'
        cmd.save(output_file)

    def color(self, tuples):
        for color, selection in tuples:
            cmd.color(color, selection)

    def color_by_element(self, selections):
        for selection in selections:
            util.cbag(selection)

    def set(self, name, value, opts=None):
        opts = opts or {}
        cmd.set(name, value, **opts)

    def reinitialize(self):
        cmd.reinitialize('everything')
        self.input_file = None

    def quit(self):
        cmd.quit()

    def run_cmds(self, commands):
        for func_name, opts in commands:
            getattr(cmd, func_name)(**opts)


class PymolColorMap:

    def __init__(self, color_map=None, default_color=None):
        self.color_map = color_map or {}
        self.default_color = default_color

    def get_color(self, key):
        if key in self.color_map:
            return self.color_map[key]
        else:
            logger.warning("Key '%s' does not exist. Default color will be used: %s." % (key, self.default_color))
            return self.default_color


def mybio_to_pymol_selection(entity):
    params = {}
    if entity.level == 'S' or entity.level == 'M':
        params[''] = 'all'
    else:
        if entity.level == 'C':
            params['chain'] = entity.id
        else:
            if entity.level == 'R':
                params['resn'] = entity.resname
                params['res'] = str(entity.id[1]) + entity.id[2].strip()
                params['chain'] = entity.get_parent().id
            else:
                if entity.level == 'A':
                    residue = entity.get_parent()
                    params['id'] = entity.serial_number
                    params['name'] = entity.name
                    params['resn'] = residue.resname
                    params['res'] = str(residue.id[1]) + residue.id[2].strip()
                    params['chain'] = residue.get_parent().id

            return (' AND ').join(['%s %s' % (k, str(v)) for k, v in params.items()])
