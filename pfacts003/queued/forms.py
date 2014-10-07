from django import forms

class EmailForm(forms.Form):
    email_to = forms.EmailField(required=False,
        label='Email to:', help_text='Email results to you (recommended: results can take an hour or longer to return)')
    email_subject = forms.CharField(required=False,
        label='Email subject:', initial='BPG Submission Results',
        help_text='Email subject line (optional)')
