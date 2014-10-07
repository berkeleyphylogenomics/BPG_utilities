import smtplib
from email.MIMEText import MIMEText
from email.MIMEMultipart import MIMEMultipart 

def send_email_from_phylofacts(recipient, text_body, html_body, subject):
    ''' Sends an email from our gmail address '''
    EMAIL_FROM = "phylofacts.webmaster@gmail.com"
    EMAIL_PASSWORD = "lanikai324c"

    text = MIMEText(text_body, 'plain')
    html = MIMEText(html_body, 'html')

    msg = MIMEMultipart('alternative')
    msg['Subject'] = subject
    msg['From'] = EMAIL_FROM
    msg['To'] = recipient
    msg.attach(text)
    msg.attach(html)

    # Sending email
    server = smtplib.SMTP('smtp.gmail.com:587')
    server.ehlo()
    server.starttls()
    server.ehlo()
    server.login(EMAIL_FROM, EMAIL_PASSWORD)
    server.sendmail(EMAIL_FROM, recipient, msg.as_string())
    server.close()
