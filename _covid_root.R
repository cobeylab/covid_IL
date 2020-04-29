#' Initializes the COVID_ROOT environment variable
covid_set_root <- function(path) {
    Sys.setenv(COVID_ROOT = normalizePath(path))
}

#' Gets a path relative to the COVID_ROOT environment variable
covid_get_path <- function(rel_path) {
    root <- Sys.getenv('COVID_ROOT')
    
    if(root == '') {
        stop(paste(
            'COVID_ROOT environment variable is not set.\n',
            'Set to repository root using covid_set_root(path).',
            '(This is because R makes it hard to find out where a file is being sourced from.'
        ))
    }
    
    file.path(root, rel_path)
}

#' Sources relative to the COVID_ROOT environment variable
covid_source <- function(rel_path) {
    source(covid_get_path(rel_path))
}
