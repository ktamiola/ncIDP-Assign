# -----------------------------------------------------------------------------
# Reposition assigned stretch of protein sequence using chemical shift
# statistics.
#
# Original Author: Tom Goddard
# Author: Kamil Tamiola
# Molecular Dynamics & NMR
# University of Groningen, The Netherlands
# k.tamiola@rug.nl
# 2nd DEC 2010 

import Tkinter
import types

import pyutil
import sequence
import shiftstats_ncIDP
import sparky
import sputil
import tkutil
import time

# -----------------------------------------------------------------------------
#
class reposition_dialog(tkutil.Dialog, tkutil.Stoppable):

  def __init__(self, session):

    self.session = session

    tkutil.Dialog.__init__(self, session.tk, 'ncIDP Repositioning Tool 1.2b/2nd DEC 2010')
    
    cm = sputil.condition_menu(session, self.top, 'Condition: ')
    cm.frame.pack(side = 'top', anchor = 'w')
    self.condition_menu = cm

    rr = tkutil.entry_row(self.top, 'Residues ',
                          ('from', '', 3), ('to', '', 3))
    self.range_variables = rr.variables
    rr.frame.pack(side = 'top', anchor = 'w')

    ol = tkutil.scrolling_list(self.top, '', 5)
    ol.frame.pack(side = 'top', fill = 'both', expand = 1)
    ol.listbox.bind('<ButtonRelease-1>', self.select_cb)
    self.offset_list = ol

    progress_label = Tkinter.Label(self.top, anchor = 'nw')
    progress_label.pack(side = 'top', anchor = 'w')

    br = tkutil.button_row(self.top,
			   ('Positions', self.positions_cb),
			   ('Move', self.move_cb),
                           ('Shifts', self.shifts_cb),
			   ('Stop', self.stop_cb),
			   ('Close', self.close_cb),
                           ('Help', sputil.help_cb(session, 'Reposition'))
			   )
    br.frame.pack(side = 'top', anchor = 'w')

    tkutil.Stoppable.__init__(self, progress_label, br.buttons[3]) 

  # ---------------------------------------------------------------------------
  # added by kT @ 2nd DEC 2010
  # this routine allows user-specified update of correction tables for ncIDP
  def update_cb(self):
    self.progress_report('Updating ncIDP...')
    time.sleep(2)
    self.progress_report('Contacting www.protein-nmr.org...')
    time.sleep(2)
    self.progress_report('Almost there...')
    #self.progress_report('Storing files in %s ' % shiftstats_ncIDP.get_sparky_home())
    

  # ---------------------------------------------------------------------------
  #
  def select_cb(self, event):

    resonance = self.offset_list.event_line_data(event)
    if type(resonance) == types.InstanceType:
      self.session.show_resonance_peak_list(resonance)

  # ---------------------------------------------------------------------------
  #
  def positions_cb(self):

    if self.read_range():
      offset_scores = self.stoppable_call(self.offset_scores, self.resonances,
                                          self.sequence, self.residue_range)
      if offset_scores != None:
        offset_scores.sort(pyutil.component_comparer(1))
        self.set_list(offset_scores, self.residue_range)

  # ---------------------------------------------------------------------------
  #
  def move_cb(self):

    offsets = self.offset_list.selected_line_data()
    if len(offsets) != 1:
      self.progress_report('Must select one list line')
      return

    offset = offsets[0]
    if type(offset) == types.IntType:
      self.reposition(self.resonances, self.sequence, offset)

  # ---------------------------------------------------------------------------
  #
  def shifts_cb(self):

    if self.read_range():
      self.offset_list.clear()
      heading = 'Group  Atom  Shift  Expected  Dev  Peaks'
      self.offset_list.heading['text'] = heading
      def nsa_value(r):
        return (r.group.number, r.group.symbol, r.atom.name)
      # here we get the resonance list
      reslist = pyutil.sort_by_function_value(self.resonances, nsa_value)
      for r in reslist:
        stats = shiftstats_ncIDP.atom_statistics(r.group.number,r.group.symbol, r.atom.name, self.RC_prediction)
        if stats:
          dev = abs(r.frequency - stats.average_shift) / stats.shift_deviation
          expected = '%6.4g %+6.1f' % (stats.average_shift, dev)
        else:
          expected = ''
        line = '%5s %5s %6.4g %14s %4d' % (r.group.name, r.atom.name,
                                           r.frequency, expected, r.peak_count)
        self.offset_list.append(line, r)
     
  # ---------------------------------------------------------------------------
  #
  def read_range(self):
    
    condition = self.condition_menu.condition()

    residue_range = (pyutil.string_to_int(self.range_variables[0].get()),
                     pyutil.string_to_int(self.range_variables[1].get()))
    if residue_range[0] == None or residue_range[1] == None:
      self.progress_report('Must select residue range.')
      return 0

    reslist = self.range_resonances(condition, residue_range)
    if len(reslist) == 0:
      self.progress_report('No residues in selected range.')
      return 0

    seq = sequence.molecule_sequence(condition.molecule)
    if seq == None:
      self.progress_report('Molecule sequence unknown.')
      return 0

    self.residue_range = residue_range
    self.resonances = reslist
    self.sequence = seq
    # Compute global ncIDP-based random coil chemical shifts
    self.RC_prediction = shiftstats_ncIDP.resonance_statistics(reslist)

    return 1
    
  # ---------------------------------------------------------------------------
  #
  def range_resonances(self, condition, residue_range):

    min, max = residue_range
    resonances = []
    for r in condition.resonance_list():
      if r.group.number >= min and r.group.number <= max:
        if r.peak_list():
          resonances.append(r)
    return resonances
    
  # ---------------------------------------------------------------------------
  #
  def set_list(self, offset_scores, residue_range):
    # Added output for clear analysis
    self.offset_list.clear()
    heading = ('Residues %d-%d\n' % residue_range +
               'Deviation Mismatches Collisions Position')
    self.offset_list.heading['text'] = heading
    #fo = open('%d-%d.dat' % residue_range  , 'w')
    for offset, mean_dev, mismatches, collisions in offset_scores:
      if mean_dev == None:
        mean_dev = 0
      location = '%d-%d' % (residue_range[0]+offset, residue_range[1]+offset)
      line = '%8.3f %8d %8d %12s' % (mean_dev, mismatches, collisions, location)
      self.offset_list.append(line, offset)
      if mean_dev != 0:
        short_line = '%8.3f  %12s\n' % (mean_dev, location)
        #fo.write(short_line)
    #fo.close() 
  # ---------------------------------------------------------------------------
  # Return list of offset, mean deviation (in SD units), mismatch triples.
  #
  def offset_scores(self, resonances, sequence, residue_range):

    first_offset = sequence.first_residue_number - residue_range[0]
    seq_length = sequence.last_residue_number-sequence.first_residue_number+1
    frag_length = residue_range[1] - residue_range[0] + 1
    offset_count = seq_length - frag_length + 1

    assigned_group_atoms = self.collision_table(resonances)
    offset_scores = []
    self.stoppable_loop('position', 1)
    for offset in range(first_offset, first_offset + offset_count):
      self.check_for_stop()
      r2ns = self.new_number_symbol(resonances, sequence, offset)
      mean_dev = self.mean_deviation(r2ns)
      mismatches = self.mismatches(r2ns)
      collisions = self.collisions(r2ns, assigned_group_atoms)
      offset_scores.append((offset, mean_dev, mismatches, collisions))
    return offset_scores

  # ---------------------------------------------------------------------------
  # Return mean deviation (in SD units)
  #
  def new_number_symbol(self, resonances, sequence, offset):

    r2ns = {}
    for r in resonances:
      num = r.group.number + offset
      if sequence.number_to_symbol.has_key(num):
        sym = sequence.number_to_symbol[num]
        r2ns[r] = (num, sym)
    return r2ns

  # ---------------------------------------------------------------------------
  # Return mean deviation (in SD units)
  #
  def mean_deviation(self, r2ns):

    total_dev = 0
    dev_count = 0
    for r, (num, sym) in r2ns.items():
      stats = shiftstats_ncIDP.atom_statistics(num,sym, r.atom.name, self.RC_prediction)
      if stats and stats.shift_deviation != 0:
        dev = abs(r.frequency - stats.average_shift) / stats.shift_deviation
        total_dev = total_dev + dev
        dev_count = dev_count + 1

    if dev_count == 0:
      return None
    return total_dev / dev_count

  # ---------------------------------------------------------------------------
  # Return number mismatches.  A mismatch is where a destination group
  # does not have an atom with the same name as the resonance being moved.
  #
  def mismatches(self, r2ns):

    mismatches = 0
    for r, (num, sym) in r2ns.items():
      atom_name = r.atom.name
      from_sym = r.group.symbol
      from_stats = shiftstats_ncIDP.atom_statistics(r.group.number,from_sym, atom_name, self.RC_prediction)
      if from_stats != None:
        stats = shiftstats_ncIDP.atom_statistics(num, sym, atom_name, self.RC_prediction)
        if stats == None:
          mismatches = mismatches + 1
    return mismatches

  # ---------------------------------------------------------------------------
  # Return table of group atoms for assigned resonances other than the
  # specified resonances but in the same condition.
  #
  def collision_table(self, resonances):
    
    group_atoms = {}

    if resonances:

      condition = resonances[0].condition
      for r in condition.resonance_list():
        if r.peak_count > 0:
          group_atoms[(r.group.name, r.atom.name)] = 1

      for r in resonances:
        ga = (r.group.name, r.atom.name)
        if group_atoms.has_key(ga):
          del group_atoms[ga]

    return group_atoms

  # ---------------------------------------------------------------------------
  # Return the number of destination resonances that are alread assigned.
  #
  def collisions(self, r2ns, assigned):

    collisions = 0
    for r, (num, sym) in r2ns.items():
      group_name = sym + repr(num) + r.group.suffix
      ga = (group_name, r.atom.name)
      if assigned.has_key(ga):
        collisions = collisions + 1
    return collisions

  # ---------------------------------------------------------------------------
  #
  def reposition(self, resonances, sequence, offset):

    new_group = {}
    for r in resonances:
      num = r.group.number + offset
      if sequence.number_to_symbol.has_key(num):
        sym = sequence.number_to_symbol[num]
      else:
        sym = 'repo'
      group_name = sym + repr(num) + r.group.suffix
      new_group[r] = group_name

    peaks = {}
    for r in resonances:
      for p in r.peak_list():
        if self.movable_peak(p, new_group):
          peaks[p] = 1
    peaks = peaks.keys()
    
    for peak in peaks:
      peak_res = peak.resonances()
      for axis in range(len(peak_res)):
        r = peak_res[axis]
        if new_group.has_key(r):
          group_name = new_group[r]
          atom_name = r.atom.name
          peak.assign(axis, group_name, atom_name)
        
  # ---------------------------------------------------------------------------
  # Is the peak assigned using resonances from table.
  #
  def movable_peak(self, peak, restable):

    for r in peak.resonances():
      if r != None and not restable.has_key(r):
        return 0
    return 1
        
  # ---------------------------------------------------------------------------
  # Start with peak and find all connected resonances and peaks.
  #
  def connected_fragment(self, peak):

    peak_table = {}
    resonance_table = {}

    self.connected_to_peak(peak, peak_table, resonance_table)

    self.peaks = peak_table.keys()
    self.resonances = resonance_table.keys()

  # ---------------------------------------------------------------------------
  # Add all connected resonances and peaks to tables.
  #
  def connected_to_peak(self, peak, peak_table, resonance_table):

    peak_table[peak] = 1
    for r in peak.resonances():
      if r and not resonance_table.has_key(r):
        self.connected_to_resonance(r, peak_table, resonance_table)

  # ---------------------------------------------------------------------------
  # Add all connected resonances and peaks to tables.
  #
  def connected_to_resonance(self, resonance, peak_table, resonance_table):

    resonance_table[resonance] = 1
    for p in resonance.peak_list():
      if not peak_table.has_key(p):
        self.connected_to_peak(p, peak_table, resonance_table)

# -----------------------------------------------------------------------------
#
def show_repositioner(session):
  sputil.the_dialog(reposition_dialog,session).show_window(1)
