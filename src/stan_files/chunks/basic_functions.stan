  vector colMeans(matrix x) {
    int ncols = cols(x);
    vector[ncols] ans;
    for (j in 1:ncols) { 
      ans[j] = mean(col(x, j));
    }
    return ans;
  }
