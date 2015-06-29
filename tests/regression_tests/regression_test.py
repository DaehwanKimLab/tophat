"""
regression_test.py

Created by Ryan Kelley
"""


import unittest
import sys
import os
import subprocess
import shutil
import filecmp
import shlex
import tempfile
#
#  Generic test case for comparing output of TopHat
#  Given a directory, it will look in that directory
#  for a file "command.txt". The test will attempt
#  to re-run the command specified in this file, and compare
#  the output to the sample output provided in the output
#  directory. The location of the "gold-standard" directory
#  is specified by the command itself. By convention, all
#  input and output files should be located beneath the parent
#  directory.
#
class TestTopHat(unittest.TestCase):
    def __init__(self, methodName, directoryName, topHatExecutable, samtoolsExecutable):
        self.directoryName = directoryName
        self.topHatExecutable = topHatExecutable
        self.samtoolsExecutable = samtoolsExecutable
        self.testOutputDirectory = tempfile.mkdtemp(dir=directoryName)
        super(TestTopHat, self).__init__(methodName=methodName)
        
    def test_output(self):
        sys.stderr.write("\nStarting test: "+self.directoryName+"\n")
	os.chdir(self.directoryName)
        #
        # Check that TopHat can be run
        # successfully
        #
        self.assertTrue(os.path.exists(self.directoryName + os.path.sep + "command.txt"), self.directoryName+": No command.txt file exists")
        commandHandle = open(self.directoryName + os.path.sep + "command.txt")
        cmd = shlex.split(commandHandle.readline())
        commandHandle.close()        
        cmd[0] = self.topHatExecutable
      
        #
        # Identify the directory that contains the output data according to the TopHat
        # command. Update the TopHat command to output the result to a new
        # temporary directory
        # 
        outputIndex = -1
        if "-o" in cmd:
            outputIndex = cmd.index("-o")
        elif "--output-dir" in cmd:
            outputIndex = cmd.index("--output-dir")
        if outputIndex > -1:
            self.goldOutputDirectory = cmd[outputIndex + 1]
            cmd = cmd[0:outputIndex] + ["--output-dir",self.testOutputDirectory] + cmd[(outputIndex + 2):]
        else:
            cmd = [cmd[0]] + ["--output-dir",self.testOutputDirectory] + cmd[1:]
            self.goldOutputDirectory = "tophat_out"
        self.fastaGenome = ""
        for prefix in cmd[-3:]:
            if os.path.exists(self.directoryName + os.path.sep + prefix + ".fa"):
                self.fastaGenome = self.directoryName + os.path.sep + prefix + ".fa"
        self.assertTrue(self.fastaGenome != "", self.directoryName+": No FASTA format genome found in genome directory")
        result = subprocess.call(cmd)    
        self.assertEqual(0,result, self.directoryName+": TopHat failed to complete successfully")

        #
        # Check that all of the appropriate output files exist
        #
        outputFiles = ["accepted_hits.bam","insertions.bed","deletions.bed", "junctions.bed"]
        for file in outputFiles:
            self.assertTrue(os.path.exists(self.testOutputDirectory + os.path.sep + file), self.directoryName+": Failed to create output file: "+file)

        #
        # Check that the BAM file can be converted to SAM
        #
        for directory in [self.testOutputDirectory, self.goldOutputDirectory]:
            cmd = [self.samtoolsExecutable]
            cmd.append("view")
            cmd.append(directory + os.path.sep + "accepted_hits.bam")
            samtoolsHandle = open(directory + os.path.sep + "accepted_hits.sam","w")
            result = subprocess.call(cmd, stdout=samtoolsHandle)
            samtoolsHandle.close()
            self.assertEqual(0, result, self.directoryName+": accepted_hits.bam failed sam/bam conversion")
	    self.assertTrue(os.path.exists(directory + os.path.sep + "accepted_hits.sam"), self.directoryName+": accepted_hits.bam failed sam/bam conversion")

        #
        # Check that samtools doesn't want any corrections
        # made in the BAM file
        #
        cmd = [self.samtoolsExecutable]
        cmd.append("calmd")
        cmd.append(self.testOutputDirectory + os.path.sep + "accepted_hits.bam")
        cmd.append(self.fastaGenome)
        nullHandle = open(os.devnull,"w")
        errorHandle = open(self.testOutputDirectory + os.path.sep + "error.txt", "w")
        result = subprocess.call(cmd, stdout = nullHandle, stderr = errorHandle)
        nullHandle.close()
        errorHandle.close()
        self.assertEqual(0, result, self.directoryName+": Unable to run \"samtools calmd\" on accepted_hits.bam")
        self.assertEqual(0, os.path.getsize(self.testOutputDirectory + os.path.sep + "error.txt"),
                         self.directoryName+": \"samtools calmd\" generated errors on accepted_hits.bam")
        

        #
        # Check that the contents of all the output files match
        # Note that the BAM file will not match because of changes
	# in the header. Instead we compare the converted SAM file
	#
        outputFiles.append("accepted_hits.sam")
        for file in outputFiles:
            result = filecmp.cmp(self.testOutputDirectory + os.path.sep + file,
                                 self.goldOutputDirectory + os.path.sep + file,
                                 shallow = False)
            self.assertTrue(result or file == "accepted_hits.bam",self.directoryName+": Validation of "+file+" failed")  

    def tearDown(self):
        if os.path.exists(self.testOutputDirectory):
            sys.stderr.write("Removing output data\n")
            shutil.rmtree(self.testOutputDirectory)            
    
def generate_TestTopHat_suite(directoryName,topHatExecutable, samtoolsExecutable):
    tests = ['test_output']
    return unittest.TestSuite([TestTopHat(methodName=test, directoryName=directoryName, topHatExecutable = topHatExecutable, samtoolsExecutable = samtoolsExecutable) for test in tests])

if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.stderr.write("Usage: python testoutput.py test_case_directory tophat_executable samtools_executable\n")
        sys.exit(-1) 
    #
    # Read and validate directory that contains
    # data for regression tests
    #
    testDirectory = sys.argv[1]
    if not os.path.isdir(testDirectory):
        sys.stderr.write(testDirectory + " should be a directory\n");
        sys.exit(-1)
    testDirectory = os.path.abspath(testDirectory)

    #
    # Read and (poorly) validate filename of
    # the TopHat executable
    #
    topHatExecutable = sys.argv[2]
    if not os.path.isfile(topHatExecutable) or (not topHatExecutable.endswith("tophat") and not topHatExecutable.endswith("tophat.py")):
        sys.stderr.write(topHatExecutable + " should be the TopHat executable\n")
        sys.exit(-1)
    topHatExecutable = os.path.abspath(topHatExecutable)

    #
    # Read and (poorly) validate filename of
    # the TopHat executable
    #
    samtoolsExecutable = sys.argv[3]
    if not os.path.isfile(samtoolsExecutable): #or not samtoolsExecutable.endswith("samtools"):
        sys.stderr.write(samtoolsExecutable + " should be the samtools executable\n")
        sys.exit(-1)
    samtoolsExecutable = os.path.abspath(samtoolsExecutable)

    testSuites = unittest.TestSuite()
    for directory in [testDirectory + os.path.sep + x for x in os.listdir(testDirectory) if os.path.isdir(testDirectory + os.path.sep + x) and x.upper().startswith("TEST")]:
        testSuites.addTests(generate_TestTopHat_suite(directory,topHatExecutable, samtoolsExecutable))
    unittest.TextTestRunner(verbosity=2).run(testSuites)
