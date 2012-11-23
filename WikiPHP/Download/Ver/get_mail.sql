select distinct (email) || ','  as dir
from xmipp_users 
where mailoption='mail' 
order by dir;

--  password usuario xmippRO 1rq6szSu
--  cat get_mail.sql |  psql -U xmippRO xmipp > list
