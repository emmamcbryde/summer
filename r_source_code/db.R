library(RSQLite)
library(XLConnect)
library(R6)



available_sheets = c('default_constants', 'country_constants', 'default_programs', 'country_programs', 'bcg_2014', 'bcg_2015',
           'bcg_2016', 'rate_birth_2014', 'rate_birth_2015', 'life_expectancy_2014', 'life_expectancy_2015',
            'notifications_2014', 'notifications_2015', 'notifications_2016', 'outcomes_2013', 'outcomes_2015',
            'mdr_2014', 'mdr_2015', 'mdr_2016', 'laboratories_2014', 'laboratories_2015', 'laboratories_2016',
            'strategy_2014', 'strategy_2015', 'strategy_2016', 'diabetes', 'gtb_2015', 'gtb_2016', 'latent_2016',
            'tb_hiv_2016', 'spending_inputs')


load_csv = function(input_path="..//xls", database_name="..//databases//Inputs.db"){
        # load csvs from input_path

        conn <- dbConnect(SQLite(), database_name)
        filenames = list.files(input_path, pattern = "\\.csv$")
        print(filenames)
        for (fname in filenames) {
          print(fname)
          tablename = sub('.csv', '', fname)
          print(tablename)
          df = read.csv(file.path(input_path, fname), header = TRUE)
          dbWriteTable(conn, tablename, df, overwrite = TRUE, row.names = FALSE)
    }
        dbDisconnect(conn)
}


load_xls = function(input_path="..//xls",  database_name="..//databases//Inputs.db" ){
        # load excel spreadsheet from input_path

        conn <- dbConnect(SQLite(), database_name)
        filenames = list.files(input_path, pattern = "\\.xlsx$")

         for (fname in filenames) {
           tablename = sub('.xlsx', '', fname)
           print(tablename)
           wb = loadWorkbook(file.path(input_path, fname))
           sheet_names = getSheets(wb)
           for (worksheet in sheet_names){
             if (worksheet %in% available_sheets){
               df = readWorksheet(wb, sheet = worksheet)
               dbWriteTable(conn, tablename, df, overwrite = TRUE, row.names = FALSE)
             }
           }
         }
        dbDisconnect(conn)
}

db_query = function(table_name, is_filter="", value="", column="*", database_name="..//databases//Inputs.db"){
        # method to query table_name

        conn <- dbConnect(SQLite(), database_name)
        query = paste("Select" , column,  "from", table_name )
        if (is_filter != "" & value != ""){
            query = paste(query, "Where", is_filter,  "=", shQuote(value)  )
        }
        res = dbFetch(dbSendQuery(conn, query))
        dbDisconnect(conn)
        return(res)
}


db_store = function(table_name, outputs, database_name="..//databases//outputs.db"){
        # method to store outputs

        conn <- dbConnect(SQLite(), database_name )
        dbWriteTable(conn, table_name, outputs, overwrite = TRUE, row.names =FALSE )
        dbDisconnect(conn)
}


#load_csv()
#load_xls()

source('sir_runner.R')
sir_model$compartment_values


db_store(table_name = 'outputs', sir_model$outputs )
res = db_query(table_name='gtb_2015', column='year')
print(res)
