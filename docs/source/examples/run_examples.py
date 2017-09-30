"""Run all examples and save their output."""

import subprocess
import os

filenames = [f for f in os.listdir('.') if os.path.isfile(f)]

for filename in filenames:
    if filename == __file__:
        # Skip this file
        continue
    if os.path.splitext(filename)[1] == '.py':
        print 'Running ' + filename
        output = subprocess.check_output(['python', filename])
        with open(os.path.splitext(filename)[0] + '_output.txt', 'w') as out_file:
            out_file.write(output)
