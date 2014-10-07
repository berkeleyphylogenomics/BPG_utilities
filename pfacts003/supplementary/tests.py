"""Supplementary Material Tests"""

from django.test import TestCase, Client

class SupplemntaryPQQTest(TestCase):
    def SetUp(self):
        self.client = Client()

    def test_url_supplementary_pqq(self):
        """Validates that /supplementary/pqq/ yields a page."""

        response = self.client.get('/supplementary/pqq/')
        self.assertEqual(response.status_code, 200)

    def test_url_supplementary_hommologs(self):
        """Validates that /supplementary/pqq/homologs_pqqa.txt yields a page."""

        response = self.client.get("/supplementary/pqq/homologs_pqqA.txt")
        self.assertEqual(response.status_code, 200)

    def test_bad_url_homologs(self):
        """Validates that /supplementary/pqq/homolgs_pqqA.txt yields a 404."""

        response = self.client.get("/supplementary/pqq/homolgs_pqqA.txt")
        self.assertEqual(response.status_code, 404)

    def test_bad_url_pqqX(self):
        """Validates that /supplementary/pqq/homologs_pqqa.txt yields a 404."""

        response = self.client.get("/supplementary/pqq/homologs_pqqa.txt")
        self.assertEqual(response.status_code, 404)

    def test_bad_url_extension(self):
        """Validates that /supplementary/pqq/homologs_pqqA.bmp yields a 404."""

        response = self.client.get("/supplementary/pqq/homologs_pqqA.bmp")
        self.assertEqual(response.status_code, 404)
