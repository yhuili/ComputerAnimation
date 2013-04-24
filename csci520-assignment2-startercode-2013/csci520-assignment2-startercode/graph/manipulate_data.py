from sys import argv

script, read_file, write_file, m = argv

print "read the data from %s" % read_file

lfemur_start = 600
lfemur_end = 800
root_start = 200
root_end = 500

in_file = open(read_file)
out_file = open(write_file, 'w')
mode = m
print "this is %s mode" % mode
line_count = 0 

while True:
    line = in_file.readline()
    if not line:
        break
    
    if (mode == 'r' and line.find('root') == 0):
    	word = line.split()
    	z = word[6]
    	line_count += 1
    	if (line_count >= root_start and line_count <= root_end):
			out_file.write(z)
			out_file.write('\n')
			print z
    elif (mode == 'l' and line.find('lfemur') == 0):
    	word = line.split()
    	x = word[1]
    	line_count += 1
    	if (line_count >= lfemur_start and line_count <= lfemur_end):
    		out_file.write(x)
			out_file.write('\n')
    		print x
	# out_file = open(write_file, 'w')

print line_count

print "all done"
out_file.close()
in_file.close()