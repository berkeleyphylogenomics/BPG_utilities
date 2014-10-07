"""Test utilities to confirm they act as expected."""

import os, subprocess
from tempfile import NamedTemporaryFile

from django.test import TestCase, Client
from django.http import Http404

from pfacts003.utils.id_patterns import is_uniprot_identifier_format,\
    is_uniprot_accession_format
from pfacts003.utils.files import serve_file
from pfacts003.utils.extract_sequence_from_fasta import extract_sequence_from_fasta
from pfacts003.utils.get_identifier_from_fasta import get_uniprot_id_from_fasta
from Bio.Alphabet import SingleLetterAlphabet
from Bio.Seq import Seq
from Bio.SeqIO import SeqRecord

class TestExtractSequenceFromFasta(TestCase):
    def test_an_input_is_required(self):
        self.assertEqual(('', '', []), extract_sequence_from_fasta(''))

    def test_an_input_cannot_simply_be_whitespace(self):
        self.assertEqual(('', '', []), extract_sequence_from_fasta('    ')) 

    def test_allowed_chars_can_be_specified(self):
        self.assertEqual(('', '', [('A', 0)]), extract_sequence_from_fasta('ABCD', allowed_chars="BCD" )) 
        self.assertEqual(('', '', [('A', 0)]), extract_sequence_from_fasta('ABCD-DC', allowed_chars="-BCD" )) 
   
    def test_the_case_of_a_sequence_does_not_matter(self):
        self.assertEqual(('', "acca", []), extract_sequence_from_fasta("acca"))

    def test_a_sequence_with_blacklist_chars(self):
         self.assertEqual(('', '',[("B", 1), ("J", 2)]), extract_sequence_from_fasta("ABJA"))

    def test_a_sequence_may_contain_newlines_tabs_and_spaces(self):
         self.assertEqual(('', "ACLACTR",[]), extract_sequence_from_fasta("ACL\nA CT\tR"))

    def test_a_sequence_may_not_contain_non_alphanumeric_chars(self):
         self.assertEqual(('', '', [('?', 3)]), extract_sequence_from_fasta("ACL?CA"))

    def test_a_sequence_may_end_with_an_asterisk(self):
         self.assertEqual(('', "ACLACTR",[]), extract_sequence_from_fasta("ACLACTR*"))

    def test_a_sequence_may_have_an_asterisk_only_at_the_end(self):
         self.assertEqual(('', '', [('*', 7)]), extract_sequence_from_fasta("ACLACTR*AC"))

    def test_a_sequence_does_not_require_a_defline(self):
         self.assertEqual(('', "ACLACTR",[]), extract_sequence_from_fasta("ACLACTR"))
         self.assertEqual(('>foo', "ACLACTR",[]), extract_sequence_from_fasta(">foo\nACLACTR"))

    def test_a_mostly_empty_defline_is_OK(self):
         self.assertEqual(('>', "ACLACTR",[]), extract_sequence_from_fasta(">\nACLACTR"))

    def test_a_sequence_cannot_just_have_a_defline(self):
         self.assertEqual(('>','', []), extract_sequence_from_fasta(">\n")) 
         self.assertEqual(('>some text here','', []), extract_sequence_from_fasta(">some text here\n"))

    def test_a_sequence_cannot_just_have_mutliple_deflines(self):
         self.assertEqual(('>defline1', '', [('>', 0)]), extract_sequence_from_fasta(">defline1\n>defline2\n"))

    def test_a_sequence_cannot_begin_with_a_gt(self):
         self.assertEqual(('>', '',[('>', 0)]), extract_sequence_from_fasta(">\n>A"))

    def test_a_defline_can_begin_with_multiple_gt_signs(self):
         self.assertEqual(('>>foo', 'ACCMNT',[]), extract_sequence_from_fasta(">>foo\nACCMNT"))

    def test_a_defline_can_contain_multiple_gt_signs(self):
         self.assertEqual(('>foo>bar', 'ACCMNT',[]), extract_sequence_from_fasta(">foo>bar\nACCMNT"))

    def test_a_sequence_may_have_at_most_one_defline(self):
         self.assertEqual(('>', '', [('>', 0), ('o', 2), ('o', 3)]), extract_sequence_from_fasta(">\n>\n\nfoo\nA"))

    def test_a_sequence_may_contain_numbers(self):
         self.assertEqual(('', 'ACCANMA',[]), extract_sequence_from_fasta("ACCANM 7 A "))

    def test_a_real_sequence_should_parse_correctly(self):
         self.assertEqual((
         '>sp|P30559|OXYR_HUMAN Oxytocin receptor OS=Homo sapiens GN=OXTR PE=1 SV=2',
         'MEGALAANWSAEAANASAAPPGAEGNRTAGPPRRNEALARVEVAVLCLILLLALSGNACVLLALRTTRQKHSRLFFFMKHLSIADLVVAVFQVLPQLLWDITFRFYGPDLLCRLVKYLQVVGMFASTYLLLLMSLDRCLAICQPLRSLRRRTDRLAVLATWLGCLVASAPQVHIFSLREVADGVFDCWAVFIQPWGPKAYITWITLAVYIVPVIVLAACYGLISFKIWQNLRLKTAAAAAAEAPEGAAAGDGGRVALARVSSVKLISKAKIRTVKMTFIIVLAFIVCWTPFFFVQMWSVWDANAPKEASAFIIVMLLASLNSCCNPWIYMLFTGHLFHELVQRFLCCSASYLKGRRLGETSASKKSNSSSFVLSHRSSSQRSCSQPSTA',[]),
         extract_sequence_from_fasta('>sp|P30559|OXYR_HUMAN Oxytocin receptor OS=Homo sapiens GN=OXTR PE=1 SV=2\nMEGALAANWSAEAANASAAPPGAEGNRTAGPPRRNEALARVEVAVLCLILLLALSGNACVLLALRTTRQKHSRLFFFMKH\nLSIADLVVAVFQVLPQLLWDITFRFYGPDLLCRLVKYLQVVGMFASTYLLLLMSLDRCLAICQPLRSLRRRTDRLAVLAT\nWLGCLVASAPQVHIFSLREVADGVFDCWAVFIQPWGPKAYITWITLAVYIVPVIVLAACYGLISFKIWQNLRLKTAAAAA\nAEAPEGAAAGDGGRVALARVSSVKLISKAKIRTVKMTFIIVLAFIVCWTPFFFVQMWSVWDANAPKEASAFIIVMLLASL\nNSCCNPWIYMLFTGHLFHELVQRFLCCSASYLKGRRLGETSASKKSNSSSFVLSHRSSSQRSCSQPSTA'))

class TestGetUniprotIdFromFasta(TestCase):
    def test_a_real_sequence_should_return_the_correct_id(self):
        self.assertEqual('sp|P30559|OXYR_HUMAN',
         get_uniprot_id_from_fasta('>sp|P30559|OXYR_HUMAN Oxytocin receptor OS=Homo sapiens GN=OXTR PE=1 SV=2\nMEGALAANWSAEAANASAAPPGAEGNRTAGPPRRNEALARVEVAVLCLILLLALSGNACVLLALRTTRQKHSRLFFFMKH\nLSIADLVVAVFQVLPQLLWDITFRFYGPDLLCRLVKYLQVVGMFASTYLLLLMSLDRCLAICQPLRSLRRRTDRLAVLAT\nWLGCLVASAPQVHIFSLREVADGVFDCWAVFIQPWGPKAYITWITLAVYIVPVIVLAACYGLISFKIWQNLRLKTAAAAA\nAEAPEGAAAGDGGRVALARVSSVKLISKAKIRTVKMTFIIVLAFIVCWTPFFFVQMWSVWDANAPKEASAFIIVMLLASL\nNSCCNPWIYMLFTGHLFHELVQRFLCCSASYLKGRRLGETSASKKSNSSSFVLSHRSSSQRSCSQPSTA')[1].id)
    def test_a_real_sequence_should_return_the_correct_description(self):
        self.assertEqual('sp|P30559|OXYR_HUMAN Oxytocin receptor OS=Homo sapiens GN=OXTR PE=1 SV=2',
         get_uniprot_id_from_fasta('>sp|P30559|OXYR_HUMAN Oxytocin receptor OS=Homo sapiens GN=OXTR PE=1 SV=2\nMEGALAANWSAEAANASAAPPGAEGNRTAGPPRRNEALARVEVAVLCLILLLALSGNACVLLALRTTRQKHSRLFFFMKH\nLSIADLVVAVFQVLPQLLWDITFRFYGPDLLCRLVKYLQVVGMFASTYLLLLMSLDRCLAICQPLRSLRRRTDRLAVLAT\nWLGCLVASAPQVHIFSLREVADGVFDCWAVFIQPWGPKAYITWITLAVYIVPVIVLAACYGLISFKIWQNLRLKTAAAAA\nAEAPEGAAAGDGGRVALARVSSVKLISKAKIRTVKMTFIIVLAFIVCWTPFFFVQMWSVWDANAPKEASAFIIVMLLASL\nNSCCNPWIYMLFTGHLFHELVQRFLCCSASYLKGRRLGETSASKKSNSSSFVLSHRSSSQRSCSQPSTA')[1].description)
    def test_a_real_sequence_should_return_the_correct_sequence(self):
        self.assertEqual('MEGALAANWSAEAANASAAPPGAEGNRTAGPPRRNEALARVEVAVLCLILLLALSGNACVLLALRTTRQKHSRLFFFMKHLSIADLVVAVFQVLPQLLWDITFRFYGPDLLCRLVKYLQVVGMFASTYLLLLMSLDRCLAICQPLRSLRRRTDRLAVLATWLGCLVASAPQVHIFSLREVADGVFDCWAVFIQPWGPKAYITWITLAVYIVPVIVLAACYGLISFKIWQNLRLKTAAAAAAEAPEGAAAGDGGRVALARVSSVKLISKAKIRTVKMTFIIVLAFIVCWTPFFFVQMWSVWDANAPKEASAFIIVMLLASLNSCCNPWIYMLFTGHLFHELVQRFLCCSASYLKGRRLGETSASKKSNSSSFVLSHRSSSQRSCSQPSTA',
         get_uniprot_id_from_fasta('>sp|P30559|OXYR_HUMAN Oxytocin receptor OS=Homo sapiens GN=OXTR PE=1 SV=2\nMEGALAANWSAEAANASAAPPGAEGNRTAGPPRRNEALARVEVAVLCLILLLALSGNACVLLALRTTRQKHSRLFFFMKH\nLSIADLVVAVFQVLPQLLWDITFRFYGPDLLCRLVKYLQVVGMFASTYLLLLMSLDRCLAICQPLRSLRRRTDRLAVLAT\nWLGCLVASAPQVHIFSLREVADGVFDCWAVFIQPWGPKAYITWITLAVYIVPVIVLAACYGLISFKIWQNLRLKTAAAAA\nAEAPEGAAAGDGGRVALARVSSVKLISKAKIRTVKMTFIIVLAFIVCWTPFFFVQMWSVWDANAPKEASAFIIVMLLASL\nNSCCNPWIYMLFTGHLFHELVQRFLCCSASYLKGRRLGETSASKKSNSSSFVLSHRSSSQRSCSQPSTA')[1].seq.tostring())

class IdentifiersAndAccessionsTest(TestCase):
    """Legacy tests for PhyloFacts (Limahuli) System"""

    def test_uniprot_accessions(self):

        """
        Verify a sample of UniProt identifiers validate
        """

        true_data = (
            'Q4JUZ8', 'O27018', 'B7GCS7', 'B5FMN5', 'A9PQS5', 'C8A7E2',
            'C5R887', 'C6QM83', 'A8XHC3', 'D1IMX6', 'D2I171', 'B3ZUA2',
            'A8FYR2', 'A8NFC2', 'Q6CHT1', 'D2QE09', 'B0HSZ7', 'B2A7G7',
            'Q4JUZ8', 'O27018', 'B7GCS7', 'B5FMN5', 'A9PQS5', 'C8A7E2',
            'C9QA98',
        )

        for data in true_data:
            self.assertTrue(is_uniprot_accession_format(data))

    def test_bogus_uniprot_accessions(self):
        """
        Verify that a sample of bad UniProt accessions fail
        """

        false_data = (
            '5R887', 'C8A72'
        )

        for data in false_data:
            self.assertFalse(is_uniprot_accession_format(data))

    def test_uniprot_identifiers(self):
        """
        Verify a sample of UniProt identifiers validate
        """

        true_data = (
            'BIOH_YERPS', 'C9M0R4_LACHE', 'A9Y305_9BACT', 'A9X730_9MAGN',
            'A4GTY6_9PROT', 'A2F8J1_TRIVA', 'INF2_MOUSE', 'RUBR1_PSEAE',
            'B2L2_SICCD', 'B32_LOXSN', 'MYG_TURTR', 'PSAI_DRIGR', 'TAL_CLOBB',
            'YBJQ_SHIB3', 'RPOC_BURCJ', 'FABZ_CAMHC',
        )

        for data in true_data:
            self.assertTrue(is_uniprot_identifier_format(data))

    def test_bogus_uniprot_identifiers(self):
        """
        Verify that a sample of bad UniProt identifiers fail
        """

        false_data = (
            'BIOH-YERPS', 'BIOH-YERP2',
        )

        for data in false_data:
            self.assertFalse(is_uniprot_identifier_format(data))

class ServeFileTest(TestCase):
    def setUp(self):
        self.disposition_string = 'attachment; filename=%s'
        self.content = "testing...1, 2, 3"
        self.filename = "test.txt"
        self.temp = NamedTemporaryFile()

        f = open(self.temp.name, 'w')
        f.write(self.content)
        f.close()

    def tearDown(self):
        #os.remove(self.temp.name)
        pass

    def test_serve_content(self):
        response = serve_file(filename=self.filename, content=self.content)
        self.assertEqual(response.content, self.content)

    def test_serve_content_as_file(self):
        disp_text = self.disposition_string % self.filename

        response = serve_file(filename=self.filename, content=self.content,
                              is_download=True)
        self.assertTrue(response.has_header('Content-Disposition'))
        self.assertEqual(response['Content-Disposition'], disp_text)

    def test_serve_content_with_custom_filename(self):
        disp_text = self.disposition_string % self.temp.name

        response = serve_file(filename=self.temp.name,
                              content=self.content,
                              is_download=True)
        self.assertTrue(response.has_header('Content-Disposition'))
        self.assertEqual(response['Content-Disposition'], disp_text)
        self.assertEqual(response.content, self.content)

    def test_serve_content_with_custom_filename(self):
        different_content = "different_content"

        response = serve_file(filename=self.temp.name,
                              content=different_content)
        self.assertEqual(response.content, different_content)

    def test_serve_file_contents(self):
        different_content = "different content"

        self.assertRaises(Http404, serve_file, 
            filename='/etc/totally_garbage_fn3434dfd3s13')

        response = serve_file(filename=self.temp.name)
        self.assertEqual(response.content, self.content)

    def test_empty_call(self):

        #Empty call serve_file()
        self.assertRaises(Http404, serve_file)
