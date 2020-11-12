
#' @export
IBTM <- function(data, k = 5, a = 50/k, b = 0.01, iter = 1000, 
                 window = 15, n_rej=10, n_part = 10, trace = FALSE,
                 biterms, detailed = FALSE, check_convergence=0, 
                 convergence_tol=1e-5, background=FALSE){
  trace <- as.integer(trace)
  n_part <- as.integer(n_part)
  n_rej <- as.integer(n_rej)
  stopifnot(k >= 1)
  stopifnot(iter >= 1)
  stopifnot(window >= 1)
  iter <- as.integer(iter)
  window <- as.integer(window)
  stopifnot(inherits(data, "data.frame"))
  if(ncol(data) == 2){
    data <- data.frame(doc_id = data[[1]], token = data[[2]], stringsAsFactors = FALSE)
  }else{
    if(!all(c("doc_id", "token") %in% colnames(data))){
      stop("please provide in data a data.frame with 2 columns as indicated in the help of BTM")
    }
  }
  data <- data[!is.na(data$doc_id) & !is.na(data$token), ]
  ## Convert tokens to integer numbers which need to be pasted into a string separated by spaces
  data$word <- factor(data$token)
  if(detailed){
    freq <- table(data$word)
    freq <- as.data.frame(freq, responseName = "freq", stringsAsFactors = FALSE)
    vocabulary <- data.frame(id = seq_along(levels(data$word)) - 1L,
                             token = levels(data$word),
                             freq = freq$freq[match(levels(data$word), freq$Var1)],
                             stringsAsFactors = FALSE)
  }else{
    vocabulary <- data.frame(id = seq_along(levels(data$word)) - 1L,
                             token = levels(data$word),
                             stringsAsFactors = FALSE)
  }
  data$word <- as.integer(data$word) - 1L

  voc <- max(data$word) + 1
  context <- split(data$word, data$doc_id)
  context <- sapply(context, FUN=function(x) paste(x, collapse = " "))

  ## Handle manual set of biterms provided by user
  if(missing(biterms)){
    biterms <- data.frame(doc_id = character(), term1 = integer(), term2 = integer(), cooc = integer(), stringsAsFactors = FALSE)
    biterms <- split(biterms, biterms$doc_id)
  }else{
    stopifnot(is.data.frame(biterms))
    if(anyNA(biterms)){
      stop("make sure there are no missing data in biterms")
    }
    if(!all(c("doc_id", "term1", "term2") %in% colnames(biterms))){
      stop("please provide in biterms a data.frame with at least 3 columns: doc_id, term1, term2, cooc - see the example in the help of BTM")
    }
    if(!all("cooc" %in% colnames(biterms))){
      biterms$cooc <- 1L
    }else{
      biterms$cooc <- as.integer(biterms$cooc)
    }
    recode <- function(x, from, to){
      to[match(x, from)]
    }
    biterms$term1 <- recode(biterms$term1, from = vocabulary$token, to = vocabulary$id)
    biterms$term2 <- recode(biterms$term2, from = vocabulary$token, to = vocabulary$id)
    if(anyNA(biterms$term1) || anyNA(biterms$term2)){
      stop("all terms in biterms should at least be available in data as well")
    }
    if(!all(biterms$doc_id %in% names(context))){
      stop("all doc_id's of the biterms should at least be available data as well")
    }
    biterms <- split(biterms, factor(biterms$doc_id, levels = names(context)), drop = FALSE)
    biterms <- lapply(biterms, FUN=function(x) as.list(x))
  }

  ## build the model
  model <- ibtm(biterms = biterms, 
                x = context, 
                K = k, 
                W = voc, 
                a = a, 
                b = b, 
                iter = iter, 
                win = window, 
                n_rej = n_rej, 
                n_part = n_part, 
                trace = as.integer(trace))
  ## make sure integer numbers are back tokens again
  rownames(model$phi) <- vocabulary$token
  ## also include vocabulary
  class(model) <- "BTM"
  if(detailed){
    model$vocabulary <- vocabulary[c("token", "freq")]
    model$biterms <- terms.BTM(model, type = "biterms")
  }
  model
}
