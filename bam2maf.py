import pysam
import argparse
from Bio.Align import MultipleSeqAlignment
from Bio.AlignIO.MafIO import MafWriter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pysam import AlignedSegment

import CigarOperation


# Adapted from https://github.com/Martinsos/edlib/issues/127
def get_pairwise_seqs(aligned_segment: AlignedSegment):
    if aligned_segment.query_name == '29637ad4-097b-4707-994d-9a4ce22c211b':
        print('stop')

    query = aligned_segment.query_sequence
    ref = aligned_segment.get_reference_sequence()
    cigar_tuples = aligned_segment.cigartuples
    ref_pos = aligned_segment.query_alignment_start
    query_pos = 0

    ref_aln = ""
    match_aln = ""
    query_aln = ""

    for operation, operation_length in cigar_tuples:
        if operation == CigarOperation.BAM_CEQUAL:
            ref_aln += ref[ref_pos: ref_pos + operation_length]
            ref_pos += operation_length
            query_aln += query[query_pos: query_pos + operation_length]
            query_pos += operation_length
            match_aln += "|" * operation_length
        elif operation == CigarOperation.BAM_CDIFF:
            ref_aln += ref[ref_pos: ref_pos + operation_length]
            ref_pos += operation_length
            query_aln += query[query_pos: query_pos + operation_length]
            query_pos += operation_length
            match_aln += "." * operation_length
        elif operation == CigarOperation.BAM_CDEL:
            ref_aln += ref[ref_pos: ref_pos + operation_length]
            ref_pos += operation_length
            query_aln += "-" * operation_length
            query_pos += 0
            match_aln += "-" * operation_length
        elif operation == CigarOperation.BAM_CINS:
            ref_aln += "-" * operation_length
            ref_pos += 0
            query_aln += query[query_pos: query_pos + operation_length]
            query_pos += operation_length
            match_aln += "-" * operation_length
        else:
            print(f"Unsupported operation {operation}")

    if len(ref_aln) != len(query_aln):
        print(aligned_segment.query_name)
        print(f'{ref_aln}\n{match_aln}\n{query_aln}')

    return ref_aln, query_aln, match_aln


# Get the number of hard clipping bases at the beginning of the alignment
def get_starting_hard_clipping_bases(aligned_segment: AlignedSegment) -> int:
    starting_hard_clipping_bases = 0
    cigar_tuples = aligned_segment.cigartuples
    if cigar_tuples is not None and len(cigar_tuples) > 0:
        first_operation, first_operation_length = cigar_tuples[0]
        if first_operation == CigarOperation.BAM_CHARD_CLIP:
            starting_hard_clipping_bases = first_operation_length
    return starting_hard_clipping_bases


# The start of the alignment is calculated as the query_alignment_start provided by pysam object + the number of
# starting hard clipped bases
def get_query_alignment_start(aligned_segment: AlignedSegment) -> int:
    starting_hard_clipping_bases = get_starting_hard_clipping_bases(aligned_segment)
    return aligned_segment.query_alignment_start + starting_hard_clipping_bases


# Receives a sam alignment and return an instance of MultipleSeqAlignment object
def create_alignment(aligned_segment: AlignedSegment) -> MultipleSeqAlignment:
    if aligned_segment.query_alignment_sequence is None:
        return None

    reference_pairwise, query_pairwise, match_pairwise = get_pairwise_seqs(aligned_segment)
    # First 's' line is the reference
    strand = 1 if aligned_segment.is_forward else -1
    reference_name = aligned_segment.reference_name
    reference_length = aligned_segment.header.get_reference_length(reference_name)
    # Size is the size of the aligning region in the source sequence.
    # This number is equal to the number of non-dash characters in the alignment text field below.
    ref_record = SeqRecord(Seq(reference_pairwise),
                           id=reference_name,
                           annotations={'start': aligned_segment.reference_start,
                                        'size':  aligned_segment.reference_length,
                                        'strand': 1,
                                        'srcSize': reference_length})

    # Second 's' line is the aligned read
    query_alignment_start = get_query_alignment_start(aligned_segment)
    seq_record = SeqRecord(Seq(query_pairwise),
                           id=aligned_segment.query_name,
                           annotations={'start': query_alignment_start,
                                        'size':  aligned_segment.query_length,
                                        'strand': strand,
                                        'srcSize': aligned_segment.infer_read_length()})
    multiple_seq_alignment = MultipleSeqAlignment([ref_record, seq_record], annotations={"score": 0})
    return multiple_seq_alignment


def bam2maf(input_file_path: str, maf_file_path: str):
    # To open a sam, the input read mode should be 'r', whereas 'rb' is needed to open a bam
    input_read_mode = 'r'
    if input_file_path.lower().endswith('.bam'):
        input_read_mode = 'rb'

    with pysam.AlignmentFile(input_file_path, input_read_mode) as input_file, open(maf_file_path, "w+") as maf_file:
        maf_writer = MafWriter(maf_file)
        maf_writer.write_header()
        for sam_entry in input_file:
            alignment = create_alignment(sam_entry)
            if alignment is not None:
                maf_writer.write_alignment(alignment)
            else:
                print(f'could not extract alignment from {sam_entry}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                    prog='bam2maf',
                    description='Transform BAM/SAM files into Multiple Alignment Format (MAF) files')

    parser.add_argument('-i', '--input', type=str, required=True, help='BAM or SAM file to be converted into MAF.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Path to store the resulting MAF file')
    args = parser.parse_args()
    bam2maf(args.input, args.output)
