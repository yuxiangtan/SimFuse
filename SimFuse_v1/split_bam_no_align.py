#   Copyright {2015} Yuxiang Tan
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

import argparse
import os
import os.path
import pysam
import sys

if pysam.__version__ == '0.7.7':
    AlignmentFile = pysam.Samfile
elif pysam.__version__ > '0.7.7':
    AlignmentFile = pysam.AlignmentFile
else:
    raise ImportError('using too old version of pysam')
    


def split_bam(filename, prefix='', output_dir='./'):
    """
    Splits filename into 4 output files, based on each read's flag:
    * Singletons
    * Secondary Alignment
    * Unmapped 
    * Paired
    """
    
    output_dir = os.path.abspath(output_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # assumes bam with 'rb'
    input_bam = AlignmentFile(os.path.abspath(filename), 'rb')
    header = input_bam.header

    output_singleton = AlignmentFile(
        os.path.join(output_dir, prefix + 'singleton.bam'), 'wb',
        header=header)
    output_paired = AlignmentFile(
        os.path.join(output_dir, prefix + 'paired.bam'), 'wb',
        header=header)
    output_unmapped = AlignmentFile(
        os.path.join(output_dir, prefix + 'unmapped.bam'), 'wb',
        header=header)
    output_secondary = AlignmentFile(
        os.path.join(output_dir, prefix + 'sec_alignment.bam'), 'wb',
        header=header)

    secondary_count = 0
    unmapped_count = 0
    singleton_count = 0
    paired_count = 0
    
    for index, read in enumerate(input_bam):
        if read.flag & 256:
            secondary_count += 1
            output_secondary.write(read)
        elif ((read.flag & 8) and (read.flag & 4)):
            unmapped_count += 1
            output_unmapped.write(read)
        elif (read.flag & 8) or (read.flag & 4):
            singleton_count += 1
            output_singleton.write(read)
        else:
            paired_count += 1
            output_paired.write(read)

    output_singleton.close()
    output_paired.close()
    output_unmapped.close()
    output_secondary.close()
    
    with open(os.path.join(output_dir, prefix + 'stats.txt'), 'w') as fout:
        fout.write('Singletons: {0:d}\n'.format(singleton_count))
        fout.write('Paired: {0:d}\n'.format(paired_count))
        fout.write('Unmapped: {0:d}\n'.format(unmapped_count))
        fout.write('Secondary Alignment: {0:d}\n'.format(secondary_count))
    print "Total number of reads: ", sum((singleton_count,
                                          paired_count,
                                          unmapped_count,
                                          secondary_count))


def main(argv):
    parser = argparse.ArgumentParser(description="""
Split bam into secondary, paired, unmapped and singleton bam files""")
    parser.add_argument('input', help='input bam file; path can be relative')
    parser.add_argument('-p', '--prefix', default='',
                        help='prefix to tack onto output file names')
    parser.add_argument('-o', '--outputdir', default='./',
                        help='which directory to write the output files to')
    opts = parser.parse_args(argv)
    split_bam(opts.input, opts.prefix, opts.outputdir)

   
if __name__ == "__main__":
    main(sys.argv[1:])
