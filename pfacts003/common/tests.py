"""Test common functionality of framework."""

from django.test import TestCase

from pfacts003.utils.testing_utils import try_path
from pfacts003.utils.id_patterns import is_uniprot_identifier_format, is_uniprot_accession_format

class LegacySystemRunningTest(TestCase):
    """Legacy tests for PhyloFacts (Limahuli) System"""

    def test_phylofacts(self):
        """
        Verify that each of these servers are running
        """
        self.assertEquals(
            try_path('phylogenomics.berkeley.edu', '/phylofacts/'), True)

class ServersRunningTest(TestCase):
    def test_makana_staging(self):
        """
        Verify staging server is running on Ohana/Makana
        """
        self.assertEquals(
            try_path('makana-test.berkeley.edu', '/'), True)

    def test_makana_production(self):
        """
        Verify production server is running on Ohana/Makana
        """
        self.assertEquals(
            try_path('makana.berkeley.edu', '/'), True)


class MiscTest(TestCase):
    """
    Verify several BPG miscellaneous pages
    """
    def test_publications(self):
        """
        Verify /publications/ link works on staging and production
        """
        self.assertEquals(
            try_path('makana.berkeley.edu', '/publications/'), True)

        self.assertEquals(
            try_path('makana-test.berkeley.edu', '/publications/'), True)

    def test_publications(self):
        """
        Verify /contact_us/ link works on staging and production
        """
        self.assertEquals(
            try_path('makana.berkeley.edu', '/contact_us/'), True)

        self.assertEquals(
            try_path('makana-test.berkeley.edu', '/contact_us/'), True)
