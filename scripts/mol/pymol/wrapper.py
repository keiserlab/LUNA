from pymol import cmd
from pymol import util
from file.util import (get_filename, get_file_format)


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

    def reset_view(self):
        cmd.delete('all')
        self.input_file = None

    def run_cmds(self, commands):
        for func_name, opts in commands:
            getattr(cmd, func_name)(**opts)
