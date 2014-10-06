Please refer to the following manual page for the manage_ohana master script.


Help on module manage_ohana:

NAME
    manage_ohana - Manage Ohana Master Script

FILE
    /clusterfs/ohana/home/glenjarvis/repositories/ohana/bpg/manage_ohana/manage_ohana.py

DESCRIPTION
    The Ohana cluster maintenance tasks, such as webserver file updates,
    database backups, etc. are encapsulated here. This serves the purpose of
    automating tasks and self-documenting procedures so we can keep more
    organized and the system in good working order.

CLASSES
    exceptions.Exception
        InvalidFunctionError
    
    class InvalidFunctionError(exceptions.Exception)
     |  Invalid function (entered as string) in menu structure
     |  
     |  You are allowed to enter a string that represents a function, instead of an
     |  actual function type, in the menu structure. However, the function should
     |  exist in globals() when the menu is attempted to be evaludated. Otherwise,
     |  this exception is raised. 
     |  
     |  Check the function name that you entered in the menu structure.
     |  
     |  Methods inherited from exceptions.Exception:
     |  
     |  __getitem__(...)
     |  
     |  __init__(...)
     |  
     |  __str__(...)

FUNCTIONS
    clear_screen()
    
    connect(...)
        connect(dsn, ...) -- Create a new database connection.
        
        This function supports two different but equivalent sets of arguments.
        A single data source name or ``dsn`` string can be used to specify the
        connection parameters, as follows::
        
            psycopg2.connect("dbname=xxx user=xxx ...")
        
        If ``dsn`` is not provided it is possible to pass the parameters as
        keyword arguments; e.g.::
        connect(dsn, ...) -- Create a new database connection.
        
        This function supports two different but equivalent sets of arguments.
        A single data source name or ``dsn`` string can be used to specify the
        connection parameters, as follows::
        
            psycopg2.connect("dbname=xxx user=xxx ...")
        
        If ``dsn`` is not provided it is possible to pass the parameters as
        keyword arguments; e.g.::
        
            psycopg2.connect(database='xxx', user='xxx', ...)
        
        The full list of available parameters is:
        
        - ``dbname`` -- database name (only in 'dsn')
        - ``database`` -- database name (only as keyword argument)
        - ``host`` -- host address (defaults to UNIX socket if not provided)
        - ``port`` -- port number (defaults to 5432 if not provided)
        - ``user`` -- user name used to authenticate
        - ``password`` -- password used to authenticate
        - ``sslmode`` -- SSL mode (see PostgreSQL documentation)
        
        - ``async`` -- if the connection should provide asynchronous API
        
        If the ``connection_factory`` keyword argument is not provided this
        function always return an instance of the `connection` class.
        Else the given sub-class of `extensions.connection` will be used to
        instantiate the connection object.
        
        :return: New database connection
        :rtype: `extensions.connection`
    
    menu(menu_data)

DATA
    BOOK_LIST_FILENAME = '/home/glenjarvis/experiments/sample_data/pfacts/...
    DB_LOG_FILEPATH = '/clusterfs/ohana/software/db//logs'
    PRODUCTION = 'production'
    SCREEN_HEIGHT = 23
    STAGING = 'staging'
    cleaning_menu = {'A': ('Clean using Apache', <function clean_webserver...
    config = <ConfigParser.SafeConfigParser instance>
    database_menu = {'B': ('Take Database Backup', <function take_db_backu...
    family_menu = {'C': ('Check Families Disk Structure', <function check_...
    main_menu = {'B': ('Book Management', 'menu', [], {'menu_data': {'C': ...
    production_webserver_menu = {'A': ('Export Everything', <function push...
    staging_webserver_menu = {'A': ('Export Everything', <function push_al...
    volume_menu = {'C': ('Check Disk Space', <function check_disk_space>, ...
    webserver_menu = {'C': ('Clean Temporary Files/Directories', 'menu', [...

