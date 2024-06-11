'''
Created on 20/01/2014

@author: MMPE

See documentation of HTCFile below

'''
import os
from collections import OrderedDict
import collections


class OrderedDict(collections.OrderedDict):
    pass

    def __str__(self):
        return "\n".join(["%-30s\t %s" % ((str(k) + ":"), str(v)) for k, v in self.items()])


def parse_next_line(lines):
    _3to2list = list(lines.pop(0).split(";"))
    line, comments, = _3to2list[:1] + [_3to2list[1:]]
    comments = ";".join(comments).rstrip()
    while lines and lines[0].lstrip().startswith(";"):
        comments += "\n%s" % lines.pop(0).rstrip()
    return line.strip(), comments


def fmt_value(v):
    try:
        if int(float(v)) == float(v):
            return int(float(v))
        return float(v)
    except ValueError:
        return v.replace("\\", "/")


class HTCContents(object):
    lines = []
    contents = None
    name_ = ""
    parent = None

    def __getitem__(self, key):
        if isinstance(key, str):
            key = key.replace(".", "/")
            if "/" in key:
                keys = key.split('/')
                val = self.contents[keys[0]]
                for k in keys[1:]:
                    val = val[k]
                return val
            return self.contents[key]
        else:
            return self.values[key]

    def __getattr__(self, *args, **kwargs):
        if args[0] in ['__members__', '__methods__']:
            # fix python2 related issue. In py2, dir(self) calls
            # __getattr__(('__members__',)), and this call must fail unhandled to work
            return object.__getattribute__(self, *args, **kwargs)
        try:
            return object.__getattribute__(self, *args, **kwargs)
        except Exception:
            k = args[0]
            if k.endswith("__1"):
                k = k[:-3]
            return self.contents[k]

    def __setattr__(self, *args, **kwargs):
        _3to2list1 = list(args)
        k, v, = _3to2list1[:1] + _3to2list1[1:]
        if k in dir(self):  # in ['section', 'filename', 'lines']:
            if isinstance(self, HTCLine) and k == 'values':
                args = k, list(v)
            return object.__setattr__(self, *args, **kwargs)
        if isinstance(v, str):
            v = [fmt_value(v) for v in v.split()]
        if not isinstance(v, HTCContents):
            if not isinstance(v, (list, tuple)):
                v = [v]
            if k in self.contents:
                self.contents[k].values = v
                return
            v = HTCLine(k, v, "")
        self.contents[k] = v
        v.parent = self

    def __delattr__(self, *args, **kwargs):
        k, = args
        if k in self:
            del self.contents[k]

    def __iter__(self):
        # mainbodies must preceed constraints
        values = ([v for v in self.contents.values() if v.name_ == 'main_body'] +
                  [v for v in self.contents.values() if v.name_ != 'main_body'])
        return iter(values)

    def __contains__(self, key):
        if self.contents is None:
            return False
        return key in self.contents

    def get(self, section, default=None):
        try:
            return self[section]
        except KeyError:
            return default

    def __call__(self, **kwargs):
        """Allow accesing one of multiple subsections with same name, e.g. the main body where name=='shaft'
        > htc.new_htc_structure.main_body(name='shaft')

        or one of multiple lines with same name, e.g. the section in c2_def where value[0]==3
        > htc.new_htc_structure.main_body.c2_def.sec(v0=3)
        """
        lst = [s for s in self.parent if s.name_ == self.name_ and (
            all([k in s and s[k][0] == v for k, v in kwargs.items()]) or
            (all([k[0] == 'v' for k in kwargs]) and all([s[int(k[1:])] == v for k, v in kwargs.items()]))
        )]
        assert len(lst) == 1
        return lst[0]

    def keys(self):
        return list(self.contents.keys())

    def _add_contents(self, contents):
        if contents.name_ not in self:
            self[contents.name_] = contents
        else:
            ending = "__2"
            while contents.name_ + ending in self:
                ending = "__%d" % (1 + float("0%s" % ending.replace("__", "")))
            self[contents.name_ + ending] = contents
        contents.parent = self

    def add_section(self, section_name, members={}, section=None, allow_duplicate=False, **kwargs):
        if isinstance(section_name, HTCSection):
            section = section_name
            section_name = section.name_
        if section_name in self and allow_duplicate is False:
            return self[section_name]
        if section_name == "output":
            section = HTCOutputSection(section_name)
        elif section_name.startswith("output_at_time"):
            section = HTCOutputAtTimeSection(section_name)
        elif section is None:
            section = HTCSection(section_name)
        self._add_contents(section)
        kwargs.update(members)
        for k, v in kwargs.items():
            section[k] = v
        return section

    def delete(self):
        keys = [k for (k, v) in self.parent.contents.items() if v == self]
        for k in keys:
            del self.parent.contents[k]

    def location(self):
        if self.parent is None:
            return os.path.basename(self.filename)
        else:
            name = [k for k in self.parent.keys() if self.parent[k] == self][0]
            return self.parent.location() + "/" + name

    def compare(self, other, compare_order=False):
        my_keys = self.keys()
        other_keys = other.keys()
        s = ""
        while my_keys or other_keys:
            if my_keys:
                if (my_keys[0] in other_keys):
                    if compare_order:
                        other_i = 0
                    else:
                        other_i = other_keys.index(my_keys[0])
                    while other_keys[other_i] != my_keys[0]:
                        s += "\n".join(["+ %s" % l for l in str(other[other_keys.pop(other_i)]
                                                                ).strip().split("\n")]) + "\n\n"

                    s += self[my_keys.pop(0)].compare(other[other_keys.pop(other_i)])

                else:
                    s += "\n".join(["- %s" % l for l in str(self[my_keys.pop(0)]).strip().split("\n")]) + "\n\n"
            else:
                s += "\n".join(["+ %s" % l for l in str(other[other_keys.pop(0)]).strip().split("\n")]) + "\n\n"
        return s


class HTCSection(HTCContents):
    end_comments = ""
    begin_comments = ""

    def __init__(self, name, begin_comments="", end_comments=""):
        self.name_ = name.strip() # strip if tabs in name somehow 
        self.begin_comments = begin_comments.strip(" \t")
        self.end_comments = end_comments.strip(" \t")
        self.contents = OrderedDict()
        self.parent = None

    @property
    def section_name(self):
        return self.name_

    @section_name.setter
    def section_name(self, value):
        self.name_ = value

    def add_line(self, name, values, comments=""):
        line = HTCLine(name, values, comments)
        self._add_contents(line)
        return line

    def __setitem__(self, key, value):
        if isinstance(value, HTCContents):
            self.contents[key] = value
            value.parent = self
        elif isinstance(value, (str, int, float)):
            self.add_line(key, [value])
        else:
            self.add_line(key, value)

    @staticmethod
    def from_lines(lines):
        line, begin_comments = parse_next_line(lines)
        name = line[6:].lower()
        if name == "output":
            section = HTCOutputSection(name, begin_comments)
        elif name.startswith("output_at_time"):
            section = HTCOutputAtTimeSection(name, begin_comments)
        else:
            section = HTCSection(name, begin_comments)
        while lines:
            if lines[0].strip() == "":
                lines.pop(0)
            if lines[0].lower().startswith("begin"):
                section._add_contents(HTCSection.from_lines(lines))
            elif lines[0].lower().startswith("end"):
                line, section.end_comments = parse_next_line(lines)
                break
            elif lines:
                section._add_contents(section.line_from_line(lines))
        else:
            raise Exception("Section '%s' has not end" % section.name_)
        return section

    def line_from_line(self, lines):
        return HTCLine.from_lines(lines)

    def __str__(self, level=0):
        s = "%sbegin %s;%s\n" % ("  " * level, self.name_, (("", "\t" + self.begin_comments)
                                                            [bool(self.begin_comments.strip())]).replace("\t\n", "\n"))
        s += "".join([c.__str__(level + 1) for c in self])
        s += "%send %s;%s\n" % ("  " * level, self.name_, (("", "\t" + self.end_comments)
                                                           [self.end_comments.strip() != ""]).replace("\t\n", "\n"))
        return s

    def get_subsection_by_name(self, name, field='name'):
        return self.get_section(name, field)

    def get_section(self, name, field='name'):
        lst = [s for s in self if field in s and s[field][0] == name]
        if len(lst) == 1:
            return lst[0]
        elif len(lst) == 0:
            raise ValueError("subsection with %s='%s' not found" % (field, name))
        else:
            raise ValueError("Multiple subsection with %s='%s' not found" % (field, name))

    def get_element(self, key, value):
        """Return subsection where subsection.<key>==value or line where line.values[key]==value"""
        if isinstance(key, int):
            lst = [s for s in self if s.values[key] == value]
        elif isinstance(key, str):
            lst = [s for s in self if key in s and s[key][0] == name]
        else:
            raise ValueError("Key argument must be int or str")
        if len(lst) == 1:
            return lst[0]
        elif len(lst) == 0:
            raise ValueError("contents with '%s=%s' not found" % (key, value))
        else:
            raise ValueError("Multiple contents with '%s=%s' not found" % (key, value))

    def copy(self):
        copy = HTCSection(name=self.name_, begin_comments=self.begin_comments, end_comments=self.end_comments)
        for k, v in self.contents.items():
            copy.contents[k] = v.copy()
        return copy


class HTCLine(HTCContents):
    values = None
    comments = ""

    def __init__(self, name, values, comments):
        if "__" in name:
            name = name[:name.index("__")]
        self.name_ = name
        self.values = list(values)
        self.comments = comments.strip(" \t")
        self.parent = None

    def __repr__(self):
        return str(self)

    def __str__(self, level=0):
        if self.name_ == "":
            return ""
        return "%s%s%s;%s\n" % ("  " * (level), self.name_,
                                ("", "\t" + self.str_values())[bool(self.values)],
                                ("", "\t" + self.comments)[bool(self.comments.strip())])

    def str_values(self):
        return " ".join([str(v) for v in self.values])

    def __getitem__(self, key):
        try:
            return self.values[key]
        except Exception:
            raise IndexError("Parameter %s does not exists for %s" % (key + 1, self.location()))

    def __setitem__(self, key, value):
        if isinstance(key, int):
            self.values[key] = value
        else:
            raise NotImplementedError

    @staticmethod
    def from_lines(lines):
        line, end_comments = parse_next_line(lines)
        if len(line.split()) > 0:
            _3to2list3 = list(line.split())
            name, values, = _3to2list3[:1] + [_3to2list3[1:]]
        else:
            name = line
            values = []

        values = [fmt_value(v) for v in values]
        return HTCLine(name, values, end_comments)

    def compare(self, other):
        s = ""
        if self.values != other.values:
            s += "\n".join(["+ %s" % l for l in str(self).strip().split("\n")]) + "\n"
            s += "\n".join(["- %s" % l for l in str(other).strip().split("\n")]) + "\n"
            s += "\n"
        return s

    def copy(self):
        return HTCLine(name=self.name_, values=self.values, comments=self.comments)


class HTCOutputSection(HTCSection):
    sensors = None

    def __init__(self, name, begin_comments="", end_comments=""):
        HTCSection.__init__(self, name, begin_comments=begin_comments, end_comments=end_comments)
        self.sensors = []

    def add_sensor(self, type, sensor, values=None, comment="", nr=None):
        values = [] if values is None else values
        self._add_sensor(HTCSensor(type, sensor, values, comment), nr)

    def _add_sensor(self, htcSensor, nr=None):
        if nr is None:
            nr = len(self.sensors)
        self.sensors.insert(nr, htcSensor)
        htcSensor.parent = self

    def line_from_line(self, lines):
        name = lines[0].split()[0].strip()

        if name in ['filename', 'data_format', 'buffer', 'time']:
            return HTCLine.from_lines(lines)
        else:
            return HTCSensor.from_lines(lines)

    def _add_contents(self, contents):
        if isinstance(contents, HTCSensor):
            self._add_sensor(contents)
        else:
            return HTCSection._add_contents(self, contents)

    def __str__(self, level=0):
        s = "%sbegin %s;%s\n" % ("  " * level, self.name_, ("", "\t" + self.begin_comments)
                                 [len(self.begin_comments.strip()) > 0])
        s += "".join([c.__str__(level + 1) for c in self])
        s += "".join([s.__str__(level + 1) for s in self.sensors])
        s += "%send %s;%s\n" % ("  " * level, self.name_, ("", "\t" + self.end_comments)
                                [self.end_comments.strip() != ""])
        return s

    def compare(self, other):
        s = HTCContents.compare(self, other)
        for s1, s2 in zip(self.sensors, other.sensors):
            s += s1.compare(s2)
        for s1 in self.sensors[len(other.sensors):]:
            s += "\n".join(["- %s" % l for l in str(s1).strip().split("\n")]) + "\n"
        for s2 in self.sensors[len(self.sensors):]:
            s += "\n".join(["- %s" % l for l in str(s2).strip().split("\n")]) + "\n"

        return s


class HTCOutputAtTimeSection(HTCOutputSection):
    type = None
    time = None

    def __init__(self, name, begin_comments="", end_comments=""):
        if len(name.split()) < 3:
            raise ValueError('"keyword" and "time" arguments required for output_at_time command:\n%s' % name)
        name, self.type, time = name.split()
        self.time = float(time)
        HTCOutputSection.__init__(self, name, begin_comments=begin_comments, end_comments=end_comments)

    def __str__(self, level=0):
        s = "%sbegin %s %s %s;%s\n" % ("  " * level, self.name_, self.type, self.time,
                                       ("", "\t" + self.begin_comments)[len(self.begin_comments.strip())])
        s += "".join([c.__str__(level + 1) for c in self])
        s += "".join([s.__str__(level + 1) for s in self.sensors])
        s += "%send %s;%s\n" % ("  " * level, self.name_, ("", "\t" + self.end_comments)
                                [self.end_comments.strip() != ""])
        return s


class HTCSensor(HTCLine):
    type = ""
    sensor = ""
    values = []

    def __init__(self, type, sensor, values, comments):
        self.type = type
        self.sensor = sensor
        self.values = list(values)
        self.comments = comments.strip(" \t")

    @staticmethod
    def from_lines(lines):
        line, comments = parse_next_line(lines)
        if len(line.split()) > 2:
            _3to2list5 = list(line.split())
            type, sensor, values, = _3to2list5[:2] + [_3to2list5[2:]]
        elif len(line.split()) == 2:
            type, sensor = line.split()
            values = []
        else:
            type, sensor, values = "", "", []

        def fmt(v):
            try:
                if int(float(v)) == float(v):
                    return int(float(v))
                return float(v)
            except ValueError:
                return v
        values = [fmt(v) for v in values]
        return HTCSensor(type, sensor, values, comments)

    def __str__(self, level=0):
        return "%s%s %s%s;%s\n" % ("  " * (level),
                                   self.type,
                                   self.sensor,
                                   ("", "\t" + self.str_values())[bool(self.values)],
                                   ("", "\t" + self.comments)[bool(self.comments.strip())])

    def delete(self):
        self.parent.sensors.remove(self)

    def compare(self, other):
        s = ""
        if self.sensor != other.sensor or self.values != other.values:
            s += "\n".join(["+ %s" % l for l in str(self).strip().split("\n")]) + "\n"
            s += "\n".join(["- %s" % l for l in str(other).strip().split("\n")]) + "\n"
            s += "\n"
        return s
