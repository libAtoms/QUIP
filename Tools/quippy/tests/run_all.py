import sys, unittest, os.path, glob


test_files = glob.glob(os.path.join(os.path.dirname(__file__), 'test*.py'))
test_mods = [ os.path.splitext(os.path.basename(t))[0] for t in test_files ]

for name in test_mods:
   __import__(name)
   globals().update(vars(sys.modules[name]))

unittest.main()
   
