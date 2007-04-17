/*--------------------------------------------------------------------------*/
struct TTimeStamp
{
	double	SecondsSinceProcessBirth;
	double	SecondsSince1904_01_01_00H00;
	struct tm
			UTCTime;
	struct tm
			LocalTime;
};

