import re

from django.http import HttpResponsePermanentRedirect
from django.conf import settings


class UrlRedirectMiddleware:
    """
    This middleware lets you match a specific url and redirect the request to a
    new url.

    You keep a tuple of url regex pattern/url redirect tuples on your site
    settings, example:

    URL_REDIRECTS = (
        (r'www\.example\.com/hello/$', 'http://hello.example.com/'),
        (r'www\.example2\.com/$', 'http://www.example.com/example2/'),
    )

    """
    def process_request(self, request):
        print "ENTRA!!!!"
        
        host = request.META['HTTP_HOST'] + request.META['PATH_INFO']
        print "HOST:", host
        
        for url_pattern, redirect_url in settings.URL_REDIRECTS:
            
            print "URL PATTERN:", url_pattern
            print "REDIRECT URL:", redirect_url 
            
            regex = re.compile(url_pattern)
            if regex.match(host):
                print "RESPONSE URL:", redirect_url
                return HttpResponsePermanentRedirect(redirect_url)