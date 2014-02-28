 CREATE TABLE xmipp_installations
 (
 id INTEGER  PRIMARY KEY,
 operating_system VARCHAR(128),
 time_installation timestamp default current_timestamp,
 remoteAddr text);

