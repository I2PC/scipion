from pyworkflow.mapper.sqlite import SqliteFlatMapper

def getSQLite(request):
    pathToSQLiteDB = request.GET.get('path')

    if pathToSQLiteDB is not None:

        if not os.path.isabs(pathToSQLiteDB):
            pathToSQLiteDB = os.path.join(os.environ['SCIPION_HOME'], pathToSQLiteDB)

        mapper = SqliteFlatMapper(pathToSQLiteDB, globals())
        records = mapper.selectAll()



        context = {
            "records": records
        }
        context = base_grid(request, context)
        context.update(csrf(request))


        return render_to_response('home/startdownload.html', context)

    else:
        redirect(download_form)

