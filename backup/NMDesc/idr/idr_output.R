
# Convert list-columns to character
plus1_control_idr2_clean <- plus1_control_idr2
plus1_control_idr2_clean[] <- lapply(plus1_control_idr2_clean, function(col) {
  if (is.list(col)) {
    sapply(col, function(x) {
      if (is.null(x)) {
        return(NA)
      } else if (length(x) == 0) {
        return("")
      } else {
        return(paste(x, collapse = ";"))
      }
    })
  } else {
    col
  }
})

# Now write to file
write.table(plus1_control_idr2_clean, file = "plus1_control_idr.txt", sep = "\t", row.names = FALSE, quote = FALSE)

