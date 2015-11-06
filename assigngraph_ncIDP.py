# -----------------------------------------------------------------------------
# Generate resonance and peak assignments for proteins using sequence and
# peak lists for several spectra.
#
# Original Author: Tom Goddard
# Author: Kamil Tamiola
# Molecular Dynamics & NMR
# University of Groningen, The Netherlands
# k.tamiola@rug.nl
# 2nd DEC 2010 

import Tkinter
import types
import axes
import expectedpeaks
import pyutil
import sequence
import shiftstats_ncIDP
import sparky
import spingraph
import sputil
import strips
import tkutil

# -----------------------------------------------------------------------------
# Specify expected peaks and propose possible resonance assignments on a
# spin graph.
#
class assigner_dialog(tkutil.Dialog, tkutil.Stoppable):

  def __init__(self, session):

    self.session = session
    self.assignments = None
    self.spin_graph_dialog = None
    self.current_resonance = None
    self.extender = None
    self.shown_group_atoms = None
    
    tkutil.Dialog.__init__(self, session.tk, 'ncIDP Spin-Graph Assigner 1.2b/2nd DEC 2010')

    lb = self.peak_list_box(self.top)
    lb.pack(side = 'top', fill = 'both', expand = 1)
    
    progress_label = Tkinter.Label(self.top, anchor = 'nw')
    progress_label.pack(side = 'top', anchor = 'w')
    br = tkutil.button_row(self.top,
			   ('Assign', self.assign_cb),
                           ('Unassign', self.unassign_cb),
                           ('Strips...', self.strip_cb),
                           ('Update', self.update_cb),
			   ('Setup...', self.setup_cb),
                           ('Stop', self.stop_cb),
                           ('Close', self.close_cb),
                           ('Help', sputil.help_cb(session, 'AssignGraph')),
			   )
    br.frame.pack(side = 'top', anchor = 'w')

    tkutil.Stoppable.__init__(self, progress_label, br.buttons[5])
# ---------------------------------------------------------------------------
  #
  def peak_list_box(self, parent):

    pl = tkutil.scrolling_list(parent, 'Possible resonance assignments', 5)
    pl.listbox.bind('<ButtonRelease-1>', self.select_peak_cb)
    pl.listbox.bind('<Double-ButtonRelease-1>', self.goto_peak_cb)
    self.peak_list = pl

    return pl.frame

  # ---------------------------------------------------------------------------
  #
  def select_peak_cb(self, event):

    ed = self.peak_list.event_line_data(event)
    if self.is_assignment(ed):
      epeak, dpeak = ed
      sputil.select_peak(dpeak)

  # ---------------------------------------------------------------------------
  #
  def goto_peak_cb(self, event):

    list_data = self.peak_list.event_line_data(event)
    if self.is_assignment(list_data):
      epeak, dpeak = list_data
      sputil.show_peak(dpeak)
    elif self.is_position(list_data):
      spectrum, freq = list_data
      pos = sputil.alias_onto_spectrum(freq, spectrum)
      sputil.show_spectrum_position(spectrum, pos)
  
  # ---------------------------------------------------------------------------
  # Make assignments associated with selected peak list line.
  # If a peak line is selected the assignment is made and any previous
  # assignments for the expected peak and the data peak are unassigned.
  # If an assignment set header line is selected then the assignments in
  # the set that are not reassignments and that have only one assignment
  # for the expected peak are made.
  #
  def assign_cb(self):

    line_data = self.peak_list.selected_line_data()
    if len(line_data) != 1:
      return

    if self.is_assignment(line_data[0]):
      self.assign_list((line_data[0],))
    elif self.is_assignment_list(line_data[0]):
      assignable = self.unique_assignments(line_data[0])
      self.assign_list(assignable)

  # ---------------------------------------------------------------------------
  # Unassign all assignments associated with selected peak list line.
  # If a peak line is selected the expected peak and data peak assignments
  # (possibly different) are unassigned.  If an assignment set header line
  # is selected then all assignments where the expected peak is assigned to
  # the data peak are unassigned.
  #
  def unassign_cb(self):

    line_data = self.peak_list.selected_line_data()
    if len(line_data) == 0 and self.current_resonance:
      r = self.current_resonance
      ed_list = self.assignments.assignments_for_resonance(r)
      self.unassign_list(ed_list)
    elif len(line_data) == 1:
      if self.is_assignment(line_data[0]):
        self.unassign_list((line_data[0],))
      elif self.is_assignment_list(line_data[0]):
        ed_list = self.exact_assignments(line_data[0])
        self.unassign_list(ed_list)
  
  # ---------------------------------------------------------------------------
  # Display strip plots strips relevant to assignment of current resonance.
  #
  def strip_cb(self):

    if self.extender == None or self.current_resonance == None:
      return

    strip_plot = strips.strip_dialog(self.session)
    epeaks = self.extender.connected_expected_peaks(self.current_resonance)
    strip_list = []
    for epeak in epeaks:
      xyz_axes = self.strip_xyz_axes(epeak)
      if xyz_axes:
        center = self.strip_center(epeak, xyz_axes)
        label = self.strip_label(epeak, xyz_axes, center)
        spectrum = epeak.peak_set.spectrum
        strip = strip_plot.spectrum_strip(spectrum, xyz_axes, center, label)
        strip.group_number = epeak.resonances[xyz_axes[0]].number
        strip_list.append(strip)

    sort_key = lambda s: (s.group_number, s.spectrum.name)
    strip_list = pyutil.sort_by_function_value(strip_list, sort_key)

    strip_plot.show_window(1)
    strip_plot.delete_all_strips()
    strip_plot.add_strips(strip_list)
  
  # ---------------------------------------------------------------------------
  #
  def strip_xyz_axes(self, epeak):

    rlist = epeak.resonances
    dim = len(rlist)
    if dim != 3:
      return None
    
    other_axes = []
    res_axis = None
    for a in range(dim):
      r = rlist[a]
      if r == self.current_resonance:
        res_axis = a
      elif self.assignments.is_resonance_assigned(r):
        other_axes.append(a)

    if len(other_axes) != 2 or res_axis == None:
      return None
    
    aX, aZ = other_axes
    aY = res_axis
    nuclei = epeak.peak_set.spectrum.nuclei
    if nuclei[aZ] == '1H' and nuclei[aX] != '1H':
      aZ, aX = other_axes
    xyz_axes = aX, aY, aZ

    return xyz_axes
  
  # ---------------------------------------------------------------------------
  #
  def strip_center(self, epeak, xyz_axes):

    aX, aY, aZ = xyz_axes
    spectrum = epeak.peak_set.spectrum
    rmin, rmax = spectrum.region
    shift = self.assignments.resonance_shift
    freq = [0,0,0]
    freq[aX] = shift(epeak.resonances[aX])
    freq[aY] = .5 * (rmin[aY] + rmax[aY])
    freq[aZ] = shift(epeak.resonances[aZ])
    center = sputil.alias_onto_spectrum(freq, spectrum)
    return center
  
  # ---------------------------------------------------------------------------
  #
  def strip_label(self, epeak, xyz_axes, center):

    aX, aY, aZ = xyz_axes
    group_name = epeak.resonances[aX].group_name
    label = '%.4g\n%.4g\n%s' % (center[aX], center[aZ], group_name)
    return label

  # ---------------------------------------------------------------------------
  #
  def is_assignment(self, list_data):

    return (type(list_data) == types.TupleType and
            len(list_data) == 2 and
            type(list_data[0]) == types.InstanceType and
            type(list_data[1]) == types.InstanceType)
  
  # ---------------------------------------------------------------------------
  #
  def is_position(self, list_data):

    return (type(list_data) == types.TupleType and
            len(list_data) == 2 and
            type(list_data[0]) == types.InstanceType and
            type(list_data[1]) == types.TupleType)
  
  # ---------------------------------------------------------------------------
  #
  def is_assignment_list(self, ed_list):

    if type(ed_list) == types.TupleType or type(ed_list) == types.ListType:
      for ed in ed_list:
        if not self.is_assignment(ed):
          return 0
    return 1
  
  # ---------------------------------------------------------------------------
  # Return assignments from list where data peak is unassigned and there is
  # only one assignment in the list for the expected peak.
  #
  def unique_assignments(self, ed_list):

    epeak_count = {}
    for epeak, dpeak in ed_list:
      if not epeak_count.has_key(epeak):
        epeak_count[epeak] = 0
      epeak_count[epeak] = epeak_count[epeak] + 1

    unique = []
    for epeak, dpeak in ed_list:
      if epeak_count[epeak] == 1 and not dpeak.is_assigned:
        unique.append((epeak, dpeak))

    return unique
  
  # ---------------------------------------------------------------------------
  # Return assignments from list where expected peak is assigned to data peak.
  #
  def exact_assignments(self, ed_list):

    exact = []
    for epeak, dpeak in ed_list:
      if self.assignments.has_assignment(epeak, dpeak):
        exact.append((epeak, dpeak))
    return exact
  
  # ---------------------------------------------------------------------------
  #
  def assign_list(self, ed_list):
      
    for epeak, dpeak in ed_list:
      self.assign_peak(epeak, dpeak)

    self.color_atoms()

  # ---------------------------------------------------------------------------
  #
  def assign_peak(self, expected_peak, data_peak):

    self.unassign_peak(expected_peak, data_peak)
    
    # Update assignments object.
    self.assignments.add_peak_assignment(expected_peak, data_peak)

    expected_peak.assign(data_peak)

    # Update spin graph to show newly assigned peak
    sgd = self.spin_graph()
    if sgd:
      sgd.add_peak(data_peak)
  
  # ---------------------------------------------------------------------------
  #
  def unassign_list(self, ed_list):

    for epeak, dpeak in ed_list:
      self.unassign_peak(epeak, dpeak)

    dpeaks = map(lambda ed: ed[1], ed_list)

    self.color_atoms()

  # ---------------------------------------------------------------------------
  #
  def unassign_peak(self, expected_peak, data_peak):

    # Update assignments object.
    self.assignments.remove_peak_assignment(expected_peak)

    # Unassign peak
    for axis in range(data_peak.spectrum.dimension):
      data_peak.assign(axis, '', '')

    sgd = self.spin_graph()
    if sgd:
      sgd.remove_peak(data_peak)
    
  # ---------------------------------------------------------------------------
  # Show list of peaks assignable with this group_atom
  #
  def atom_select_cb(self, group_atom):

    asys = self.assignments.assignment_system
    if asys.sequence == None:
      return

    resonance = asys.sequence.group_atom_to_resonance[group_atom]
    self.show_assignments(resonance)

  # ---------------------------------------------------------------------------
  # Show sets of possible assignments for this resonance.
  #
  def show_assignments(self, resonance):

    self.extender = assignment_extender(self.assignments)
    extensions = self.extender.assignment_extensions(resonance, self)
    extensions = pyutil.sort_by_function_value(extensions,
                                               self.extension_sort_key)
    extensions.reverse()

    self.peak_list.clear()
    self.peak_list.heading['text'] = self.assignment_list_heading(resonance)
    for extension in extensions:
      self.show_assignment_extension(extension)

    self.current_resonance = resonance
    
  # ---------------------------------------------------------------------------
  #
  def assignment_list_heading(self, resonance):

    heading = 'Resonance %s %s' % (resonance.group_name, resonance.atom_name)
    if self.assignments.is_resonance_assigned(resonance):
      shift = self.assignments.resonance_shift(resonance)
      dev_text = self.deviation_text(resonance, shift)
      heading = heading + ' %6.4g %s' % (shift, dev_text)
    return heading
    
  # ---------------------------------------------------------------------------
  #
  def extension_sort_key(self, extension):

    ed_list = extension.peak_assignments
    r_shifts = extension.resonance_shifts

    epeak_table = {}
    for epeak, dpeak in ed_list:
      epeak_table[epeak] = 1

    no_reassignment_epeak_table = {}
    for epeak, dpeak in ed_list:
      if (not dpeak.is_assigned or
          self.assignments.has_assignment(epeak, dpeak)):
        no_reassignment_epeak_table[epeak] = 1

    return (len(epeak_table), len(no_reassignment_epeak_table))
  
  # ---------------------------------------------------------------------------
  #
  def show_assignment_extension(self, extension):

    rshifts = extension.resonance_shifts
    rtext = self.resonances_text(rshifts)
    alist = extension.peak_assignments
    self.peak_list.append(rtext, alist)

    def assignment_sort(ed):
      return (ed[0].peak_set.spectrum.name, ed[0].assignment_text())
    ualist = map(lambda epeak: (epeak, None),
                 extension.unassigned_expected_peaks())
    ed_list = pyutil.sort_by_function_value(alist + ualist, assignment_sort)
                 
    prev_epeak = None
    for epeak, dpeak in ed_list:
      show_assignment = (epeak != prev_epeak)
      freq = extension.expected_peak_frequency(epeak)
      atext = self.assignment_text(epeak, dpeak, freq, show_assignment)
      if dpeak == None:
        list_data = (epeak.peak_set.spectrum, freq)
      else:
        list_data = (epeak, dpeak)
      self.peak_list.append(atext, list_data)
      prev_epeak = epeak

    self.peak_list.append('', None)
    
  # ---------------------------------------------------------------------------
  #
  def resonances_text(self, rshifts):

    rtext = ''
    for res, shift in rshifts.items():
      if rtext:
        rtext = rtext + '   '
      rtext = rtext + self.resonance_shift_text(res, shift)
    return rtext
  
  # ---------------------------------------------------------------------------
  #
  def resonance_shift_text(self, resonance, shift):

    gname = resonance.group_name
    aname = resonance.atom_name
    dev = self.deviation_text(resonance, shift)
    return '%s %s %6.3f %s computed using ncIDP' % (gname, aname, shift, dev)
  
  # ---------------------------------------------------------------------------
  #
  def deviation_text(self, resonance, shift):

    # Added by kT on 30th Nov 2010 follwing the suggestions by Tom Goddard
    # from our e-mail correspondence (16th Sep 2010)
    asys = self.assignments.assignment_system
    self.RC_prediction = shiftstats_ncIDP.sequence_resonance_statistics(asys.sequence)
    # The end of the addition
    stats = shiftstats_ncIDP.atom_statistics(resonance.number, resonance.symbol, resonance.atom_name, self.RC_prediction )
    if stats == None or stats.shift_deviation == 0:
      return ''

    dev = (shift - stats.average_shift) / stats.shift_deviation
    return '(%+.3f SD)' % dev

  # ---------------------------------------------------------------------------
  #
  def assignment_text(self, epeak, dpeak, freq, show_assignment):

    sp = epeak.peak_set.spectrum
    if show_assignment:
      spectrum = sp.name
      assignment = epeak.assignment_text()
    else:
      spectrum = ''
      assignment = ''

    if dpeak:
      position = self.position_text(freq, dpeak.frequency)
      height = sputil.peak_height(dpeak) / sp.noise
      assigned = dpeak.assignment
    else:
      position = ''
      pos = sputil.alias_onto_spectrum(freq, sp)
      height = sp.data_height(pos) / sp.noise
      assigned = ''
      
    text = ('%10s %20s %21s %5.1f %s'
            % (spectrum, assignment, position, height, assigned))
    
    return text
      
  # ---------------------------------------------------------------------------
  #
  def position_text(self, efreq, dfreq):

    ptext = ''
    for a in range(len(dfreq)):
      offset = dfreq[a] - efreq[a]
      if offset >= 0:
        component = ' +%5.3f' % offset
      else:
        component = ' %6.3f' % offset
      ptext = ptext + component

    return ptext

  # ---------------------------------------------------------------------------
  #
  def setup_cb(self):

    setup_dialog = sputil.the_dialog(assigner_setup_dialog, self.session)
    if self.assignments:
      asys = self.assignments.assignment_system
    else:
      asys = None
    setup_dialog.set_parent_dialog(self, asys, self.update_assignment_system)
    setup_dialog.show_window(1)

  # ---------------------------------------------------------------------------
  #
  def update_cb(self):

    if self.assignments:
      asys = self.assignments.assignment_system
      asys.reread_data_peaks()
      self.update_assignment_system(asys)

  # ---------------------------------------------------------------------------
  #
  def update_assignments(self, assignment_system):

    a = assignments(assignment_system)
    a.use_spectrum_assignments(self)
    self.assignments = a

  # ---------------------------------------------------------------------------
  #
  def update_graph_peaks(self):

    sgd = self.spin_graph()
    if sgd:
      sgd.unplot_all_spectra()
      for asp in self.assignments.assignment_system.assignable_spectra:
        s = asp.spectrum
        self.progress_report('Showing peaks for ' + s.name)
        sgd.plot_spectrum(s)
      self.color_atoms()

  # ---------------------------------------------------------------------------
  #
  def update_assignment_system(self, assignment_system):

    self.update_assignments(assignment_system)
    self.update_graph_atoms()
    self.update_graph_peaks()

    setup_dialog = sputil.the_dialog(assigner_setup_dialog, self.session)
    setup_dialog.show_settings(assignment_system)

    self.progress_report('')

  # ---------------------------------------------------------------------------
  #
  def spin_graph(self):
    
    sgd = self.spin_graph_dialog
    if sgd == None or sgd.is_window_destroyed():
      sgd = spingraph.spin_graph_dialog(self.session)
      sgd.pointer_action.add_atom_select_callback(self.atom_select_cb)
      self.spin_graph_dialog = sgd
      self.shown_group_atoms = None
    return self.spin_graph_dialog

  # ---------------------------------------------------------------------------
  #
  def update_graph_atoms(self):
    
    sgd = self.spin_graph()
    if sgd == None:
      return

    sgd.show_window(1)

    asys = self.assignments.assignment_system
    if asys.sequence == None:
      return

    galist = pyutil.attribute_values(asys.assignable_resonances, 'group_atom')
    if galist != self.shown_group_atoms:
      self.shown_group_atoms = galist
      self.progress_report('Showing spin graph atoms')
      groups_per_row = asys.sequence.last_residue_number + 1
      sgd.set_graph(galist, groups_per_row, 0, None)
      sgd.graph.update_what_to_show(show_line_label = 0,
                                    lines_to_groups = 0,
                                    shade_by_range = 0)
      sgd.graph.show_atoms(1)
      sgd.graph.full_view('y')

  # ---------------------------------------------------------------------------
  #
  def color_atoms(self):

    sgd = self.spin_graph()
    if sgd == None:
      return
    
    self.progress_report('Coloring atoms')
    for r in self.assignments.assigned_resonances():
      sgd.shade_group_atom(r.group_atom, '')
    for r in self.assignments.unassigned_resonances():
      sgd.shade_group_atom(r.group_atom, 'gray25')
    assignable = self.assignments.assignable_resonances()
    for r in assignable:
      sgd.shade_group_atom(r.group_atom, 'gray50')
    self.progress_report('%d assignable resonances' % len(assignable))

# -----------------------------------------------------------------------------
# Select spectra, expected peaks, and tolerances for assigner_dialog.
#
class assigner_setup_dialog(tkutil.Settings_Dialog, tkutil.Stoppable):

  def __init__(self, session):

    self.session = session

    tkutil.Settings_Dialog.__init__(self, session.tk, 'Assignment Graph Setup')

    sc = self.spectrum_choice_table(self.top)
    sc.pack(side = 'top', anchor = 'w')

    te = self.tolerance_entry(self.top)
    te.pack(side = 'top', anchor = 'w')
    
    progress_label = Tkinter.Label(self.top, anchor = 'nw')
    progress_label.pack(side = 'top', anchor = 'w')
    
    br = tkutil.button_row(self.top,
                           ('Ok', self.ok_cb),
			   ('Apply', self.apply_cb),
                           ('Stop', self.stop_cb),
                           ('Close', self.close_cb),
                           ('Help', sputil.help_cb(session, 'AssignGraph')),
			   )
    br.frame.pack(side = 'top', anchor = 'w')

    tkutil.Stoppable.__init__(self, progress_label, br.buttons[2])

  # ---------------------------------------------------------------------------
  #
  def spectrum_choice_table(self, parent):

    headings = ('Use spectra  ', 'Expected peaks  ', 'Pattern axes')
    st = sputil.spectrum_table(self.session, parent, headings,
                               self.add_spectrum, self.remove_spectrum)
    st.spectrum_to_checkbutton = {}
    st.chosen_spectra = []
    st.spectrum_epeak_menus = {}
    st.axis_order_menu = {}
    self.spectrum_table = st
    
    spectra = self.session.project.spectrum_list()
    for spectrum in spectra:
      st.add_spectrum(spectrum)

    return st.frame

  # ---------------------------------------------------------------------------
  #
  def add_spectrum(self, spectrum, table, row):

    pat_name, pat_axes = expectedpeaks.recall_pattern(spectrum)

    #
    # Make spectrum checkbutton
    #
    cb = tkutil.checkbutton(table.frame, spectrum.name, 0)
    cb.button['selectcolor'] = sputil.spectrum_color(spectrum)
    choose_cb = pyutil.precompose(sputil.choose_spectrum_cb, spectrum,
                                  table.chosen_spectra)
    cb.add_callback(choose_cb)
    cb.button.grid(row = row, column = 0, sticky = 'w')
    if pat_name:
      cb.set_state(1)
    table.spectrum_to_checkbutton[spectrum] = cb

    #
    # Make peak pattern menu
    #
    epeak_types = expectedpeaks.expected_peak_descriptions.keys()
    pm = tkutil.option_menu(table.frame, '', epeak_types)
    pm.frame.grid(row = row, column = 1, sticky = 'w')
    table.spectrum_epeak_menus[spectrum] = pm

    #
    # Get default spectrum peak pattern and axis order
    #
    if pat_name == None:
      pat_name = expectedpeaks.default_spectrum_pattern_name(spectrum)
      if pat_name == None:
        pat_name = pm.get()
      pat_list = expectedpeaks.expected_peak_descriptions[pat_name]
      pat_nuclei = expectedpeaks.pattern_nuclei(pat_list[0])
      pat_axes = pyutil.order(pat_nuclei, spectrum.nuclei)

    #
    # Make axis order menu
    #
    aom = axes.axis_order_menu(table.frame, spectrum.nuclei,
                               initial_order = pat_axes)
    aom.frame.grid(row = row, column = 2, sticky = 'w')
    table.axis_order_menu[spectrum] = aom

    #
    # Register a callback to limit the menu choices for axis order.
    #
    def restrict_axis_order_cb(pat_name, aom=aom):
      import expectedpeaks
      pat_list = expectedpeaks.expected_peak_descriptions[pat_name]
      pat_nuclei = expectedpeaks.pattern_nuclei(pat_list[0])
      aom.restrict_menu_permutations(pat_nuclei)
    pm.add_callback(restrict_axis_order_cb)

    pm.set(pat_name)
    
  # ---------------------------------------------------------------------------
  #
  def remove_spectrum(self, spectrum, table):

    cb = table.spectrum_to_checkbutton[spectrum]
    cb.set_state(0)
    cb.button.destroy()
    del table.spectrum_to_checkbutton[spectrum]

    table.spectrum_epeak_menus[spectrum].frame.destroy()
    del table.spectrum_epeak_menus[spectrum]

    table.axis_order_menu[spectrum].frame.destroy()
    del table.axis_order_menu[spectrum]
    
  # ---------------------------------------------------------------------------
  #
  def tolerance_entry(self, parent):

    t = tkutil.entry_row(parent, 'Tolerances (ppm)',
                         (' 1H', '.02', 5),
                         (' 13C', '.4', 5),
                         (' 15N', '.2', 5))
    self.tolerances = {'1H': t.variables[0],
                       '13C': t.variables[1],
                       '15N': t.variables[2]}

    return t.frame

  # ---------------------------------------------------------------------------
  #
  def read_assignable_spectra(self):

    tolerance_table = self.read_tolerances()
    
    asplist = []
    st = self.spectrum_table
    spectra = st.chosen_spectra
    for spectrum in spectra:
      tolerances = pyutil.table_values(tolerance_table, spectrum.nuclei)
      pat_name = st.spectrum_epeak_menus[spectrum].get()
      pat_axes = st.axis_order_menu[spectrum].axis_order()
      asp = self.stoppable_call(expectedpeaks.assignable_spectrum_from_pattern,
                                spectrum, tolerances, pat_name, pat_axes, self)
      if asp == None:
        break
      asplist.append(asp)

    return asplist

  # ---------------------------------------------------------------------------
  #
  def read_tolerances(self):

    ttable = {}
    for nucleus, tol in self.tolerances.items():
      ttable[nucleus] = pyutil.string_to_float(tol.get(), 0)
    return ttable

  # ---------------------------------------------------------------------------
  # Update checkbuttons to reflect currently selected spectra
  #
  def show_settings(self, asys):

    st = self.spectrum_table
    s2cb = st.spectrum_to_checkbutton
    for cb in s2cb.values():
      cb.set_state(0)

    if asys:
      s2epm = st.spectrum_epeak_menus
      s2aom = st.axis_order_menu
      for asp in asys.assignable_spectra:
        s = asp.spectrum
        if s2cb.has_key(s):
          s2cb[s].set_state(1)
          if hasattr(asp, 'expected_peak_pattern_name'):
            s2epm[s].set(asp.expected_peak_pattern_name)
          if hasattr(asp, 'pattern_axis_order'):
            s2aom[s].set_axis_order(asp.pattern_axis_order)
          for axis in range(s.dimension):
            var = self.tolerances[s.nuclei[axis]]
            cur_tol = pyutil.string_to_float(var.get(), 0)
            new_tol = asp.tolerances[axis]          
            if new_tol != cur_tol:
              var.set('%.3g' % new_tol)
          
  # ---------------------------------------------------------------------------
  #
  def get_settings(self):

    asp_list = self.read_assignable_spectra()
    asys = expectedpeaks.assignment_system(asp_list)
    return asys

# -----------------------------------------------------------------------------
# Record peak and resonance assignments.
# This duplicates Sparky's ability to manage peak and resonance assignments.
# It is needed so that sets of conflicting assignments can be remembered.
#
class assignments:

  def __init__(self, assignment_system):

    self.assignment_system = assignment_system
    self.expected_peak_to_data_peak = {}

    #
    # Resonance assignment information is derived from peak assignments.
    # These tables are maintained for efficiency.
    #
    self.resonance_to_assignments = {}
    self.resonance_to_shift = {}
    
  # ---------------------------------------------------------------------------
  #
  def add_peak_assignment(self, expected_peak, data_peak):

    self.remove_peak_assignment(expected_peak)
    self.expected_peak_to_data_peak[expected_peak] = data_peak
    for r in expected_peak.resonances:
      self.add_resonance_assignment(r, expected_peak, data_peak)
    
  # ---------------------------------------------------------------------------
  #
  def remove_peak_assignment(self, expected_peak):

    if self.expected_peak_to_data_peak.has_key(expected_peak):

      data_peak = self.expected_peak_to_data_peak[expected_peak]
      del self.expected_peak_to_data_peak[expected_peak]
      for r in expected_peak.resonances:
        self.remove_resonance_assignment(r, expected_peak)
        
  # ---------------------------------------------------------------------------
  #
  def add_resonance_assignment(self, resonance, expected_peak, data_peak):

    if not self.resonance_to_assignments.has_key(resonance):
      self.resonance_to_assignments[resonance] = {}
    self.resonance_to_assignments[resonance][expected_peak] = data_peak
    self.resonance_to_shift[resonance] = None
        
  # ---------------------------------------------------------------------------
  #
  def remove_resonance_assignment(self, resonance, expected_peak):

    if self.resonance_to_assignments.has_key(resonance):
      e2d = self.resonance_to_assignments[resonance]
      if e2d.has_key(expected_peak):
        del e2d[expected_peak]
        self.resonance_to_shift[resonance] = None
        
  # ---------------------------------------------------------------------------
  #
  def resonance_shift(self, resonance):

    r2s = self.resonance_to_shift
    if not r2s.has_key(resonance) or r2s[resonance] == None:
      r2s[resonance] = self.compute_shift(resonance)
    return r2s[resonance]
        
  # ---------------------------------------------------------------------------
  #
  def compute_shift(self, resonance):

    if not self.resonance_to_assignments.has_key(resonance):
      return None

    total_shift = 0
    shift_count = 0
    for epeak, dpeak in self.resonance_to_assignments[resonance].items():
      for axis in epeak.resonance_axes(resonance):
        total_shift = total_shift + dpeak.frequency[axis]
        shift_count = shift_count + 1

    if shift_count == 0:
      return None

    return total_shift / shift_count
    
  # ---------------------------------------------------------------------------
  #
  def has_assignment(self, expected_peak, data_peak):

    return (self.expected_peak_to_data_peak.has_key(expected_peak) and
            self.expected_peak_to_data_peak[expected_peak] == data_peak)
    
  # ---------------------------------------------------------------------------
  #
  def is_peak_assigned(self, expected_peak):

    return self.expected_peak_to_data_peak.has_key(expected_peak)
    
  # ---------------------------------------------------------------------------
  #
  def is_resonance_assigned(self, resonance):

    r2a = self.resonance_to_assignments
    return (r2a.has_key(resonance) and len(r2a[resonance]) > 0)
    
  # ---------------------------------------------------------------------------
  #
  def is_resonance_unassigned(self, resonance):

    r2a = self.resonance_to_assignments
    return (not r2a.has_key(resonance) or len(r2a[resonance]) == 0)
    
  # ---------------------------------------------------------------------------
  #
  def any_resonances_assigned(self, resonances):

    r2a = self.resonance_to_assignments
    for r in resonances:
      if r2a.has_key(r) and r2a[r]:
        return 1
    return 0
    
  # ---------------------------------------------------------------------------
  #
  def use_spectrum_assignments(self, stoppable):

    for asp in self.assignment_system.assignable_spectra:
      stoppable.progress_report('Calculating possible assignments for ' +
                                asp.spectrum.name)
      for dpeak in asp.data_peaks:
        if dpeak.is_assigned:
          epeak = asp.expected_peak_for_assigned_data_peak(dpeak)
          if epeak:
            self.add_peak_assignment(epeak, dpeak)
    
  # ---------------------------------------------------------------------------
  # Return list of assignments (ie epeak, dpeak pairs) for peaks involving
  # the specified resonance.
  #
  def assignments_for_resonance(self, resonance):

    r2a = self.resonance_to_assignments
    if r2a.has_key(resonance):
      return r2a[resonance].items()
    return []
    
  # ---------------------------------------------------------------------------
  #
  def assigned_resonances(self):

    assigned = []
    for r, e2d in self.resonance_to_assignments.items():
      if e2d:
        assigned.append(r)
    return assigned
    
  # ---------------------------------------------------------------------------
  #
  def unassigned_resonances(self):

    assignable = self.assignment_system.assignable_resonances
    assigned = self.assigned_resonances()
    return pyutil.subtract_lists(assignable, assigned)
    
  # ---------------------------------------------------------------------------
  #
  def peak_assignments(self):

    return self.expected_peak_to_data_peak.items()
    
  # ---------------------------------------------------------------------------
  #
  def assigned_peak_count(self):

    return len(self.expected_peak_to_data_peak)

  # ---------------------------------------------------------------------------
  #
  def existing_assignments(self, epeaks):

    ed_list = []
    ep2dp = self.expected_peak_to_data_peak
    for epeak in epeaks:
      if ep2dp.has_key(epeak):
        dpeak = ep2dp[epeak]
        ed_list.append((epeak, dpeak))
    return ed_list
  
  # ---------------------------------------------------------------------------
  #
  def unassigned_peak_resonances(self, epeak):

    return filter(self.is_resonance_unassigned, epeak.resonances)
    
  # ---------------------------------------------------------------------------
  #
  def assignable_resonances(self):

    assigned = {}
    for r in self.assigned_resonances():
      assigned[r] = 1

    r2ep = self.assignment_system.resonance_to_expected_peaks
    assignable = {}
    for resonance in assigned.keys():
      for epeak in r2ep[resonance]:
        for r in epeak.resonances:
          if not assigned.has_key(r):
            assignable[r] = 1
    
    return assignable.keys()

# -----------------------------------------------------------------------------
# Find possible assignment extensions for assigning a given resonance.
#
class assignment_extender:

  def __init__(self, assignments):

    self.assignments = assignments

  # ---------------------------------------------------------------------------
  #
  def assignment_extensions(self, resonance, stoppable):

    a = self.assignments
    if a.is_resonance_assigned(resonance):
      elist = [self.peak_assignment_extension(resonance, stoppable)] 
      return elist

    elist = []
    epeaks = self.connected_expected_peaks(resonance)
    print 'Connected expected peaks ', len(epeaks)
    epeaks = self.minimally_unassigned_expected_peaks(epeaks)
    print 'Min unassigned expected peaks ', len(epeaks)
    cache = {}
    for epeak in epeaks:
      ur = tuple(a.unassigned_peak_resonances(epeak))
      if not cache.has_key(ur):
        cache[ur] = self.assignable_expected_peaks(ur)
      assignable = cache[ur]
      freq = map(a.resonance_shift, epeak.resonances)
      dpeaks = epeak.matching_data_peaks(freq, stoppable)
      print 'Matching data peaks ', len(dpeaks)
      for dpeak in dpeaks:
        rshifts = self.unassigned_shift_table(epeak, dpeak)
        e = assignment_extension(assignable, rshifts, a, stoppable)
        elist.append(e)

    #
    # Eliminate extensions with exactly the same peak assignments
    #
    etable = {}
    for e in elist:
      alist = list(e.peak_assignments)
      alist.sort()
      etable[tuple(alist)] = e
    elist = etable.values()

    return elist
  
  # ---------------------------------------------------------------------------
  # Find additional peak assignments for resonance that is already assigned.
  #
  def peak_assignment_extension(self, resonance, stoppable):
    
    a = self.assignments
    epeaks = self.assignable_expected_peaks([resonance])
    rshifts = {}
    e = assignment_extension(epeaks, rshifts, a, stoppable)
    e.include_existing_assignments()
    return e

  # ---------------------------------------------------------------------------
  #
  def minimally_unassigned_expected_peaks(self, epeaks):

    if len(epeaks) == 0:
      return []
    
    a = self.assignments
    urtable = pyutil.table_map(a.unassigned_peak_resonances, epeaks)
    urcounts = map(lambda seq: len(seq), urtable.values())
    min_unassigned = min(urcounts)
    eplist = []
    for epeak, urlist in urtable.items():
      if len(urlist) == min_unassigned:
        eplist.append(epeak)
    return eplist

  # ---------------------------------------------------------------------------
  #
  def connected_expected_peaks(self, resonance):

    connected = []
    a = self.assignments
    asys = a.assignment_system
    epeaks = asys.resonance_to_expected_peaks[resonance]
    for epeak in epeaks:
      if a.any_resonances_assigned(epeak.resonances):
        connected.append(epeak)
    return connected

  # ---------------------------------------------------------------------------
  #
  def assignable_expected_peaks(self, resonances):

    rtable = {}
    for r in resonances:
      rtable[r] = 1

    asys = self.assignments.assignment_system
    assignable = {}
    for r in resonances:
      epeaks = asys.resonance_to_expected_peaks[r]
      for epeak in epeaks:
        if not assignable.has_key(epeak):
          if self.are_resonances_assigned(epeak.resonances, rtable):
            assignable[epeak] = 1

    return assignable.keys()
    
  # ---------------------------------------------------------------------------
  #
  def are_resonances_assigned(self, resonances, rtable):

    a = self.assignments
    for r in resonances:
      if not (rtable.has_key(r) or a.is_resonance_assigned(r)):
        return 0
    return 1

  # ---------------------------------------------------------------------------
  #
  def unassigned_shift_table(self, epeak, dpeak):

    a = self.assignments
    rshifts = {}
    for axis in range(len(epeak.resonances)):
      r = epeak.resonances[axis]
      if a.is_resonance_unassigned(r):
        rshifts[r] = dpeak.frequency[axis]
    return rshifts

# -----------------------------------------------------------------------------
# Find a list of peak assignments consistent with a set of resonance
# assignments that extend already accepted assignments.
#
class assignment_extension:
    
  def __init__(self, epeaks, rshifts, assignments, stoppable):

    self.resonance_shifts = rshifts
    self.expected_peaks = epeaks
    self.assignments = assignments
    self.peak_assignments = \
      self.consistent_assignments(epeaks, stoppable)

  # ---------------------------------------------------------------------------
  #
  def consistent_assignments(self, epeaks, stoppable):

    ed_list = []
    for epeak in epeaks:
      dpeaks = self.consistent_data_peaks(epeak, stoppable)
      for dpeak in dpeaks:
        ed_list.append((epeak, dpeak))
    return ed_list
  
  # ---------------------------------------------------------------------------
  #
  def consistent_data_peaks(self, epeak, stoppable):

    freq = self.expected_peak_frequency(epeak)
    dpeaks = epeak.matching_data_peaks(freq, stoppable)
    return dpeaks

  # ---------------------------------------------------------------------------
  #
  def expected_peak_frequency(self, epeak):

    rshifts = self.resonance_shifts
    assignments = self.assignments
    freq = []
    for r in epeak.resonances:
      if rshifts.has_key(r):
        freq.append(rshifts[r])
      else:
        freq.append(assignments.resonance_shift(r))
    return tuple(freq)
  
  # ---------------------------------------------------------------------------
  #
  def unassigned_expected_peaks(self):

    assigned = map(lambda ed: ed[0], self.peak_assignments)
    all = self.expected_peaks
    unassigned = pyutil.subtract_lists(all, assigned)
    return unassigned
  
  # ---------------------------------------------------------------------------
  #
  def include_existing_assignments(self):

    a = self.assignments
    existing_assignments = a.existing_assignments(self.expected_peaks)
    alist = self.peak_assignments + existing_assignments
    self.peak_assignments = pyutil.unique_elements(alist)
    
# -----------------------------------------------------------------------------
#
def show_assignment_system(session, asys):
  d = sputil.the_dialog(assigner_dialog, session)
  d.show_window(1)
  d.update_assignment_system(asys)

# -----------------------------------------------------------------------------
#
def show_assigner(session):
  sputil.the_dialog(assigner_dialog,session).show_window(1)
