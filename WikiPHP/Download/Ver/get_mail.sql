select distinct (email) || ','  as dir
from xmipp_users 
where mailoption='mail' 
order by dir;

--  password usuario UUUUUUU PPPPPPPP
--  cat get_mail.sql |  psql -U UUUUUUU xmipp > list
