#!/usr/bin/python

from django.test import TestCase, Client

from pfacts003.phylofacts import forms

class TestNewSequenceSearchForm(TestCase):
    def test_an_empty_sequence_generates_the_correct_response(self):
        f = forms.NewSequenceSearchForm('')
        self.assertEqual("This field is required" in str(f), True)
       
    def test_a_nucleotide_sequence_generates_the_correct_response(self):
        f = forms.NewSequenceSearchForm({'input' : 'acgt'})
        self.assertEqual("Only amino acid sequences are accepted" in str(f), True)

    def test_illegal_characters_are_reported_in_the_error_message(self):
        f = forms.NewSequenceSearchForm({'input' : 'abcgt'})
        self.assertEqual("Error: illegal characters were found at: position 2,(&#39;b&#39;)" in str(f), True)

    def test_a_submission_that_consists_only_of_a_defline_generates_the_correct_response(self):
        f = forms.NewSequenceSearchForm({'input' : '>some defline\n'})
        self.assertEqual("You did not submit a sequence" in str(f), True)
