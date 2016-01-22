from django.core.mail import send_mail


def _send_mail(subject, message, sender, to, fail_silently):
    # Wrapper to DJANGO email function
    # We might want to test if email config is set or something else.

    # If to is not a list
    to = formatTo(to)

    send_mail(subject, message, sender, to, fail_silently)


def formatTo(to):

    if not isinstance(to, list):
        to = [to]

    return to


def subscribeToUsersList(emailToSubscribe):
    from pyworkflow.web.pages.settings import EMAIL_CONF
    message = 'Please, subscribe me to the scipion users list. Sent from Scipion web site on behalf of '\
              + emailToSubscribe
    to = EMAIL_CONF['EMAIL_FOR_SUBSCRIPTION']

    subject = 'Please, subscribe ' + emailToSubscribe
    sender = emailToSubscribe

    _send_mail(subject, message, sender, to, True)


def validateEmail(email):
    from django.core.validators import validate_email
    from django.core.exceptions import ValidationError
    try:
        validate_email(email)
        return True
    except ValidationError:
        return False