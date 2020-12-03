  vector colMeans(matrix x) {
    int ncols = cols(x);
    vector[ncols] ans;
    for (j in 1:ncols) { 
      ans[j] = mean(col(x, j));
    }
    return ans;
  }

  vector bbootColMeans_rng(matrix x) {
    int nrow = rows(x);
    int ncol = cols(x);
    vector[nrow] weight_vec;
    vector[ncol] ans;
    weight_vec = dirichlet_rng(rep_vector(1, nrow)); # ./ nrow;
    for (j in 1:ncol) {
      ans[j] = dot_product(weight_vec, col(x, j));
    }
    return ans;
  }

  matrix uncollapse_matrix(matrix x, int[] w) {
    matrix[sum(w), cols(x)] big_x;
    int i = 1;
    for (n in 1:rows(x)) {
      for (w_i in 1:w[n]) {
        // fill rows of big_x with row n a total of w[n] times
        big_x[i,] = x[n, ];
        i += 1;
      }
    }
    return big_x;
  }
