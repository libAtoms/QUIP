#!/usr/bin/env python

import sys, os, optparse, re

p = optparse.OptionParser(usage="%prog infile outfile [force] [virial] [other options]")

p.add_option("-n","--no_old_density", action="store_true", default=False, help="""If true, don't try to use old density""")
(opts, args) = p.parse_args()

if len(args) < 2:
   p.error("At least infile and outfile are required")

if not os.environ.has_key("PROFESS"):
   p.error("PROFESS environment variable with path of executable is required")
profess=os.environ['PROFESS']
if not os.environ.has_key("PROFESS_TEMPLATE"):
   profess_template='profess'
else:
   profess_template=os.environ['PROFESS_TEMPLATE']

infile=args[0]
outfile=args[1]

do_force=False
do_virial=False
energy_field="energy"
for i in range(2,len(args)):
   m=re.match("force=?(\S+)?", args[i])
   if m is not None:
      do_force=True
      if m.group(1) is not None:
	 force_field=m.group(1)
      else:
	 force_field="force"
   m=re.match("virial=?(\S+)?", args[i])
   if m is not None:
      do_virial=True
      if m.group(1) is not None:
	 virial_field=m.group(1)
      else:
	 virial_field="virial"
   m=re.match("energy=?(\S+)?", args[i])
   if m is not None:
      if m.group(1) is not None:
	 energy_field=m.group(1)
      else:
	 energy_field="energy"

if do_virial:
   sys.stderr.write("No support for virial yet\n")
   sys.exit(1)

# make and cd run dir
if not os.path.exists("profess_driver_run"):
   os.mkdir("profess_driver_run")

os.chdir("profess_driver_run")

# create input from template
f_inpt_template = open("../%s.inpt.template" % profess_template)
f_inpt_out = open("profess_driver.inpt","w")
l_print=[]
for l in f_inpt_template:
   print_this_line=True
   if re.match("\s*geometryfile", l) or \
      re.match("\s*calc\s+for", l) or \
      re.match("\s*calc\s+str", l) or \
      re.match("\s*tran", l) or \
      re.match("\s*rhof", l):
      print_this_line=False
   if print_this_line:
      f_inpt_out.write(l)
   else:
      f_inpt_out.write("## "+l)
f_inpt_template.close()

f_inpt_out.write("geometryfile profess_driver.ion\n")
if (do_force):
   f_inpt_out.write("calc for\n")
if (do_virial):
   f_inpt_out.write("calc str\n")
f_inpt_out.write("tran on\n")
f_inpt_out.write("prin den\n")
if not opts.no_old_density and os.path.exists("profess_driver.0.den"):
      f_inpt_out.write("rhof profess_driver.0.den")
f_inpt_out.close()

# create geometry from template
f_ion_out = open("profess_driver.ion","w")
f_input_pos = open("../"+args[0])
pos_lines=f_input_pos.readlines()
N_at=int(pos_lines[0])
if not re.search("Properties=species:S:1:pos:R:3", pos_lines[1]):
   sys.stderr.write("Need Properties=species:S:1:pos:R:3\n")
   sys.exit(3)
m=re.search('Lattice=["\'{]\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*["\'}]',pos_lines[1])
if m is None:
   sys.stderr.write("Failed to parse Lattice from %s\n" % args[0])
   sys.exit(2)
lattice_str=m.group(0)
f_ion_out.write("%BLOCK LATTICE_CART\n")
f_ion_out.write("%s %s %s\n" % (m.group(1), m.group(2), m.group(3)))
f_ion_out.write("%s %s %s\n" % (m.group(4), m.group(5), m.group(6)))
f_ion_out.write("%s %s %s\n" % (m.group(7), m.group(8), m.group(9)))
f_ion_out.write("%END BLOCK LATTICE_CART\n")
f_ion_out.write("%BLOCK POSITIONS_CART\n")
for i in range(2,len(pos_lines)):
   fields=pos_lines[i].split()
   f_ion_out.write ("%s %s %s %s\n" % (fields[0], fields[1], fields[2], fields[3]))
   pos_lines[i]="%s %s %s %s" % (fields[0], fields[1], fields[2], fields[3])
f_ion_out.write("%END BLOCK POSITIONS_CART\n")
f_ion_template = open("../%s.ion.template" % profess_template)
just_got_species_pot=False
in_species_pot=False
in_no_print=False
for l in f_ion_template:
   if re.match("\s*%BLOCK\s+(LATTICE|POSITIONS)",l):
	 in_no_print=True
   if re.match("\s*%END\s+BLOCK\s+(LATTICE|POSITIONS)",l):
	 in_no_print=False
   if re.match("\s*%BLOCK\s+SPECIES_POT",l):
      just_got_species_pot=True
   if re.match("\s*%END\s+BLOCK\s+SPECIES_POT",l):
      in_species_pot=False
   if (not in_no_print):
      f_ion_out.write(l)
   if (just_got_species_pot):
      just_got_species_pot=False
      in_species_pot=True
   elif (in_species_pot):
      fields=l.split()
      if not re.search('/', fields[1].strip()):
	 os.system("cp ../%s %s" % (fields[1].strip(), fields[1].strip()))
f_ion_template.close()
f_ion_out.close()

os.system(profess+" profess_driver")

f_profess_out=open("profess_driver.trans")
profess_out_ls=f_profess_out.readlines()
energy=float(profess_out_ls[0].strip())
force_str=[]
if (do_force):
   for i in range(N_at):
      force_str.append(profess_out_ls[2+i].strip())

f_pos_output=open("../"+args[1],"w")
f_pos_output.write("%d\n" % N_at)
f_pos_output.write("%s %s=%f Properties=species:S:1:pos:R:3" % (lattice_str, energy_field, energy))
if (do_force):
   f_pos_output.write(":%s:R:3" % force_field)
if (do_virial):
   f_pos_output.write(" %s=\"%f %f %f   %f %f %f   %f %f %f\"" % virial_field)
f_pos_output.write("\n")
for i in range(N_at):
   f_pos_output.write(pos_lines[2+i])
   if (do_force):
      f_pos_output.write(" "+force_str[i])
   f_pos_output.write("\n")
f_pos_output.close()
