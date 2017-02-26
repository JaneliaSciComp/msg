
"""
Simple script to send emails from the command line.
Example call:
python msg/send_email.py -e 10.11.5.23 -t test@janelia.hhmi.org -s test -b test
"""

import os
import sys
import optparse
import datetime

try:
    #http://pypi.python.org/pypi/mailer
    from mailer import Mailer
    from mailer import Message
except ImportError:
    #This isn't a critical service so die without throwing an error
    print "Email alerts are not available"
    sys.exit(0)

def send_email(email_host, to, cc, subject, body):
    """send the PDF as an email """   
    def tolist(email_str):
        email_str = email_str or ''
        email_str = email_str.replace(',',';')
        if ';' in email_str:
            return email_str.split(';')
        else:
            return [email_str]
    message = Message(From=tolist(to)[0], To=tolist(to), CC=tolist(cc), charset="utf-8")
    message.Subject = subject
    #message.Html = """This email uses <strong>HTML</strong>!"""
    message.Body = body
    #message.attach(filename=report_path, cid="Scheduled_Report.pdf")
    sender = Mailer(email_host)
    sender.send(message)
    
def main():
    """Parse command line args, and call appropriate functions."""
    usage="""\
usage: %prog [options]
"""
    parser = optparse.OptionParser(usage=usage)
    #Other option types are int and float, string is default.
    #Note there is also a default parameter.
    parser.add_option('-e','--host',dest="host",type="string")
    parser.add_option('-t','--to',dest="to",type="string")
    parser.add_option('-c','--cc',dest="cc",type="string")
    parser.add_option('-s','--subject',dest="subject",default="",type="string")
    parser.add_option('-b','--body',dest="body",default="",type="string")
    opts,args=parser.parse_args() #Args taken from sys.argv[1:] by default, parsed using GNU/POSIX syntax.
    if not opts.host:
        parser.error("An SMTP host address is required")

    #Update template slots in subject
    # Subject can contain these string which will be replaced by current value:
    #{DATE}, {LANG}, {DATETIME}
    opts.subject = opts.subject.replace('{DATE}', datetime.date.today().isoformat())
    opts.subject = opts.subject.replace('{DATETIME}', datetime.datetime.today().isoformat())
    opts.subject = opts.subject.replace('{LANG}', 'Python')

    send_email(opts.host, opts.to, opts.cc, opts.subject, opts.body)
    
if __name__=='__main__':
    main()