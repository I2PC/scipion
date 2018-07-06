function [OUT1 OUT2 OUT3] = projUnitBall(IN1,IN2,IN3)

if nargin == 1
    OUT1 = IN1 ./ max( 1,sqrt(IN1.^2) );
end
%
if nargin == 3
    OUT1 = IN1 ./ max( 1,sqrt(IN1.^2 + IN2.^2 + IN3.^2) );
    OUT2 = IN2 ./ max( 1,sqrt(IN1.^2 + IN2.^2 + IN3.^2) );
    OUT3 = IN3 ./ max( 1,sqrt(IN1.^2 + IN2.^2 + IN3.^2) );
end

end