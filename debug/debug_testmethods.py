import testmethods as tm

test = tm.NsabaTestMethods()
test.aba_load('../data_dir')
test.ns_load('../data_dir')

test.cross_check_coords()

