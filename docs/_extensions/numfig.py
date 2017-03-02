##
## @file    numfig.py
##
## @brief   For numbering figures in Sphinx
##
## @version $Id: numfig.py 10 2012-10-07 14:39:47Z jrobcary $
##
## Started from https://bitbucket.org/arjones6/sphinx-numfig/wiki/Home
##
## Copyright &copy; 2005-2012, Tech-X Corporation, Boulder, CO
## Free for any use whatsoever.
##

from docutils.nodes \
  import figure, caption, Text, reference, raw, SkipNode, Element
from sphinx.roles import XRefRole

#
# Element classes
#
class page_ref(reference):
  pass

class num_ref(reference):
  pass

#
# Visit/depart functions
#
# Why is SkipNode raised?
#
def skip_page_ref(self, node):
  raise SkipNode

def skip_num_ref(self, node):
  raise SkipNode

def latex_visit_page_ref(self, node):
  self.body.append("\\pageref{%s:%s}" % (node['refdoc'], node['reftarget']))
  raise SkipNode

def latex_visit_num_ref(self, node):
  fields = node['reftarget'].split('#')
  if len(fields) > 1:
    label, target = fields
    ref_link = '%s:%s' % (node['refdoc'], target)
    latex = "\\hyperref[%s]{%s \\ref*{%s}}" % (ref_link, label, ref_link)
    self.body.append(latex)
  else:
    self.body.append('\\ref{%s:%s}' % (node['refdoc'], fields[0]))
  raise SkipNode

def latex_depart_num_ref(self, node):
  pass

def html_visit_num_ref(self, node):
  fields = node['reftarget'].split('#')
  if len(fields) > 1:
    label, target = fields
    target_file = ''
    if node['refdoc']==target_file:
# Target file and curent file are the same
      link = "%s.html#%s" %(node['refdoc'], target.lower())
    else:
      link = "%s.html#%s" %(target_file, target.lower())
    html = '<a href="%s">%s</a>' %(link,  label)
    self.body.append(html)
  else:
    self.body.append('<a href="%s.html">%s</a>' % (node['refdoc'], fields[0]))

def html_depart_num_ref(self, node):
  pass

def compute_numfig_fignums(app, doctree):
# Generate figure numbers for each figure
  env = app.builder.env
  i = getattr(env, 'i', 1)
  figids = getattr(env, 'figids', {})
  figid_docname_map = getattr(env, 'figid_docname_map', {})
  for figure_info in doctree.traverse(figure):
    if app.builder.name != 'latex' and app.config.numfig_number_figures:
      for cap in figure_info.traverse(caption):
        cap[0] = Text("%s %d: %s" % \
          (app.config.numfig_figure_caption_prefix, i, cap[0]))
    for id in figure_info['ids']:
      figids[id] = i
      figid_docname_map[id] = env.docname
    i += 1
  env.figid_docname_map = figid_docname_map
  env.i = i
  env.figids = figids

def insert_numfig_links(app, doctree, docname):
# Replace numfig nodes with links
  figids = app.builder.env.figids
  if app.builder.name != 'latex':
    for ref_info in doctree.traverse(num_ref):

      if '#' in ref_info['reftarget']:
        label, target = ref_info['reftarget'].split('#')
        labelfmt = label + " %d"
      else:
        labelfmt = '%d'
        target = ref_info['reftarget']

      if target not in figids:
        continue

      if app.builder.name == 'html':
        target_doc = app.builder.env.figid_docname_map[target]
        link = "%s#%s" % (app.builder.get_relative_uri(docname, target_doc),
                 target)
        html = '<a href="%s">%s</a>' % (link, labelfmt %(figids[target]))
        ref_info.replace_self(raw(html, html, format='html'))
      else:
        ref_info.replace_self(Text(labelfmt % (figids[target])))

def setup(app):

# Are these used?
  app.add_config_value('numfig_number_figures', True, True)
  app.add_config_value('numfig_figure_caption_prefix', "Figure", True)

  app.add_node(page_ref,
    text=(skip_page_ref, None),
    html=(skip_page_ref, None),
    latex=(latex_visit_page_ref, None))

  app.add_role('page', XRefRole(nodeclass=page_ref))

  app.add_node(num_ref,
    text=(skip_num_ref, None),
    html=(html_visit_num_ref, html_depart_num_ref),
    latex=(latex_visit_num_ref, latex_depart_num_ref))

  app.add_role('num', XRefRole(nodeclass=num_ref))

  app.connect('doctree-read', compute_numfig_fignums)
  app.connect('doctree-resolved', insert_numfig_links)

