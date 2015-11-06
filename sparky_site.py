# -----------------------------------------------------------------------------
# Sparky site initialization file.  The initialize_session routine in this
# file is invoked when a Sparky session is started.  We add menu entries
# for standard extensions.
#
import sys
import sparky

# -----------------------------------------------------------------------------
#
def initialize_session(session):

  sme = standard_menu_entries()
  add_menu_entries(session, sme)

# -----------------------------------------------------------------------------
#
def standard_menu_entries():

  assignment_menu = (
    ('ad', 'Assignment distances',      ('distance','show_dialog')),
    ('aa', 'AutoAssign',                ('autoassign','show_dialog')),
    ('cs', 'Chemical shift plot',       ('chemshift','show_shifts_dialog')),
    ('md', 'Mirror peak dialog',	('mirror','show_mirror_dialog')),
    ('mp', 'Mirror peak show',          ('mirror','show_mirror_peak')),
    ('na', 'Noesy assignment',          ('noesyassign','show_dialog')),
    ('rs', 'Reposition sequence',       ('reposition','show_repositioner')),
    ('sg', 'Spin graph',		('spingraph','show_spin_graph')),
    ('ga', 'Spin graph assigner',       ('assigngraph','show_assigner')),
    )

  integration_menu = (
    ('cl', 'Copy linewidths',           ('copylinewidth',
                                         'copy_linewidths_and_positions')),
    ('lp', 'Linewidth plot',            ('linewidthplot',
                                         'show_linewidth_dialog')),
    ('rh', 'Relaxation peak heights',   ('relax','show_peak_heights')),
    ('ve', 'Volume errors',             ('volumeerror',
                                         'show_volume_error_dialog')),
    )

  molecule_menu = (
    ('ax', 'Atom name translation',     ('atoms','show_translations_dialog')),
    ('km', 'Chimera molecule view',     ('chimeraview','show_chimera_dialog')),
    ('mc', 'Midas constraints',         ('midasconstraint',
                                         'show_constraint_dialog')),
    ('ma', 'Midas atom picking',	('midaspick','show_atom_pick_dialog')),
    ('sq', 'Molecule sequence',         ('sequence','show_sequence_dialog')),
    ('pn', 'PDB atom names',            ('pdb','show_pdb_atom_dialog')),
    )

  peak_menu = (
    ('hc', 'HC peaks',                  ('hcpeaks','show_dialog')),
    ('SL', 'List selected peaks',  	('peaklist','show_selected_peaks')),
    ('mf', 'MARDIGRAS format',          ('mardigras','show_dialog')),
    ('LT', 'Peak list',         	('peaklist','show_spectrum_peaks')),
    ('pb', 'Peak table',		('peaktable','show_peak_table')),
    ('rp', 'Read peak list',       	('readpeaks','read_peak_list')),
    ('kr', 'Restricted peak pick',      ('restrictedpick','show_dialog')),
    ('mv', 'Shift resonances',          ('movepeaks','show_move_peak_dialog')),
    ('xe', 'XEASY, DYANA format',       ('xeasy','show_dialog')),
    ('xf', 'XPLOR format',              ('xplor','show_dialog')),
    )

  spectrum_menu = (
    ('al', 'Align spectrum',            ('align','show_shift_dialog')),
    ('cx', 'CORMA spectrum',            ('cormaspectrum',
                                         'show_corma_spectrum')),
    ('fm', 'Open multiple files', 	('openspectra','show_file_dialog')),
    ('la', 'Spectrum labelled axis',    ('axes','show_attached_axis_dialog')),
    ('rm', 'Spectrum region RMSD',      ('regionrmsd',
                                         'show_region_rmsd_dialog')),
    )

  fold_menu = (
    ('f1', 'Add w1 sweepwidth',         ('foldspectrum',
                                         'fold_selected_spectrum', 0, 1)),
    ('f2', 'Add w2 sweepwidth',         ('foldspectrum',
                                         'fold_selected_spectrum', 1, 1)),
    ('f3', 'Add w3 sweepwidth',         ('foldspectrum',
                                         'fold_selected_spectrum', 2, 1)),
    ('f4', 'Add w4 sweepwidth',         ('foldspectrum',
                                         'fold_selected_spectrum', 3, 1)),
    ('F1', 'Subtract w1 sweepwidth',    ('foldspectrum',
                                         'fold_selected_spectrum', 0, -1)),
    ('F2', 'Subtract w2 sweepwidth',    ('foldspectrum',
                                         'fold_selected_spectrum', 1, -1)),
    ('F3', 'Subtract w3 sweepwidth',    ('foldspectrum',
                                         'fold_selected_spectrum', 2, -1)),
    ('F4', 'Subtract w4 sweepwidth',    ('foldspectrum',
                                         'fold_selected_spectrum', 3, -1)),
    )

  view_menu = (
    ('cv', 'Center view setup',         ('centerview','center_view_setup')),
    ('sp', 'Strip plot',                ('strips','show_strip_plot')),
    )

  top_menu = (
    ('py', 'Python shell',              ('pythonshell','show_python_shell')),
    )

  # Added by kT
  rug_nmr_menu = (
    ('RS', 'ncIDP Repositioner 1.2b',        ('reposition_ncIDP', 'show_repositioner')),
    ('SG', 'ncIDP Spin Graph 1.2b',                ('assigngraph_ncIDP','show_assigner')),
    )

  #
  # The menu list consists of pairs having a menu path and entry list.
  # The menu path describes where under the Extension menu the entries should
  # appear.  An empty menu path puts the entries in the Extension menu.  If
  # the menu path is 'Spectrum' the entries are put in a cascaded menu called
  # Spectrum.  A menu path of 'Spectrum/Fold' puts the entries in a Fold menu
  # cascaded under the spectrum menu.
  # 
  menus = (
    ('',                        top_menu),
    ('Assignment',              assignment_menu),
    ('Integration',             integration_menu),
    ('Molecule',                molecule_menu),
    ('Peak',                    peak_menu),
    ('Spectrum',                spectrum_menu),
    ('Spectrum/Fold spectrum',  fold_menu),
    ('View',                    view_menu),
    ('ncIDP',                   rug_nmr_menu),
    )

  return menus

# -----------------------------------------------------------------------------
#
def add_menu_entries(session, menu_list):

  for menu_path, entries in menu_list:
    for accelerator, menu_text, function_spec in entries:

      if menu_path:
        path = menu_path + '/' + menu_text
      else:
        path = menu_text

      def func(module_name = function_spec[0],
               function_name = function_spec[1],
               session = session,
               extra_args = function_spec[2:]):
        import pythonshell
        pythonshell.invoke_module_function(session, module_name, function_name,
                                           (session,) + extra_args)

      session.add_command(accelerator, path, func)
