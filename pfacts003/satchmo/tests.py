"""
Test common functionality of framework.
"""

import httplib

from django.test import TestCase

from pfacts003.utils.testing_utils import try_path


class SatchmoRunningTest(TestCase):

    def test_production(self):
        """
        Verify that each of these servers are running
        """
        self.assertEquals(try_path('makana.berkeley.edu', '/satchmo'), True)
        self.assertEquals(try_path('makana.berkeley.edu', '/q/satchmo'), True)

    def test_staging(self):
        """
        Verify that each of these servers are running
        """
        self.assertEquals(
            try_path('makana-test.berkeley.edu', '/satchmo'), True)
        self.assertEquals(
            try_path('makana-test.berkeley.edu', '/q/satchmo'), True)

    def test_help_about_pages(self):
        """
        Verify About pages don't give errors
        """
        self.assertEquals(
            try_path('makana-test.berkeley.edu', '/satchmo/about/'), True)
        self.assertEquals(
            try_path('makana.berkeley.edu', '/satchmo/about/'), True)
        self.assertEquals(
            try_path('makana-test.berkeley.edu', '/satchmo/help/'), True)
        self.assertEquals(
            try_path('makana.berkeley.edu', '/satchmo/help/'), True)
