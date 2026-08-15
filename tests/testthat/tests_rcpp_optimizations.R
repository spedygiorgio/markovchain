context("C++ Backend Optimizations and Bugfixes")

test_that("generatorToTransitionMatrix correctly handles absorbing states (Issue #211)", {
  # Creiamo una matrice generatrice con uno stato assorbente ("C")
  gen_mat <- matrix(c(
    -2,  2,  0,
     1, -1,  0,
     0,  0,  0 
  ), nrow = 3, byrow = TRUE)
  rownames(gen_mat) <- colnames(gen_mat) <- c("A", "B", "C")
  
  # Chiamiamo la funzione (se hai usato il punto e creato il wrapper R usa la versione senza punto)
  trans_mat <- generatorToTransitionMatrix(gen_mat, byrow = TRUE)
  
  # 1. Non ci devono essere NaN derivati da 0/0
  expect_false(any(is.nan(trans_mat)))
  
  # 2. La probabilità di transizione per gli stati normali deve essere corretta (-gen/gen)
  expect_equal(trans_mat["A", "B"], 1)
  
  # 3. Lo stato assorbente "C" deve avere 1 sulla diagonale e 0 altrove
  expect_equal(trans_mat["C", "C"], 1)
  expect_equal(trans_mat["C", "A"], 0)
  expect_equal(trans_mat["C", "B"], 0)
})

test_that("clean_nas parses memory efficiently without dropping valid sequence data", {
  # clean_nas viene chiamata da createSequenceMatrix per pulire i missing values
  seq_with_nas <- c("a", "b", NA, "a", "c", "c", NA, "b")
  
  # sanitize = FALSE per mantenere i conteggi grezzi
  mat <- createSequenceMatrix(seq_with_nas, sanitize = FALSE)
  
  # Verifica che le dimensioni della matrice considerino solo "a", "b", "c"
  expect_equal(nrow(mat), 3)
  expect_equal(colnames(mat), c("a", "b", "c"))
  
  # "NA" non deve essere diventato uno stato della catena
  expect_false("NA" %in% rownames(mat))
  
  # Verifica che le transizioni valide attorno ai NA siano state estratte correttamente
  expect_equal(mat["a", "b"], 1)
  expect_equal(mat["a", "c"], 1)
  expect_equal(mat["c", "c"], 1)
  
  # Transizioni inesistenti
  expect_equal(mat["b", "a"], 0)
})

test_that("sortByDimNames (hash map optimization) sorts matrices preserving data integrity", {
  # sortByDimNames è invocata da inferHyperparam per riordinare la matrice di input.
  # Creiamo una matrice volutamente non in ordine alfabetico.
  states_unsorted <- c("Z", "A", "M")
  P_unsorted <- matrix(c(
    0.1, 0.8, 0.1,  # Riga Z
    0.2, 0.6, 0.2,  # Riga A
    0.3, 0.3, 0.4   # Riga M
  ), nrow = 3, byrow = TRUE, dimnames = list(states_unsorted, states_unsorted))
  
  scale_vec <- c(10, 10, 10)
  
  # inferHyperparam scala la matrice e la ordina alfabeticamente usando sortByDimNames
  result <- inferHyperparam(transMatr = P_unsorted, scale = scale_vec)
  sorted_mat <- result$scaledInference
  
  # 1. Verifica che righe e colonne siano state riordinate alfabeticamente
  expected_names <- c("A", "M", "Z")
  expect_equal(rownames(sorted_mat), expected_names)
  expect_equal(colnames(sorted_mat), expected_names)
  
  # 2. Verifica che i valori della matrice abbiano seguito correttamente il riordino O(1)
  # In P_unsorted, P["Z", "A"] era 0.8. Moltiplicato per lo scale(10) deve essere 8.
  expect_equal(sorted_mat["Z", "A"], 8)
  
  # P["A", "M"] era 0.2 -> Scalato: 2
  expect_equal(sorted_mat["A", "M"], 2)
  
  # P["M", "Z"] era 0.4 -> Scalato: 4
  expect_equal(sorted_mat["M", "Z"], 4)
})