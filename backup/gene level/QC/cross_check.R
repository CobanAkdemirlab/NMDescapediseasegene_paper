#fs2
aenmd_NMD = data.frame(NMDesc_can = fs2@elementMetadata@listData[["res_aenmd"]]@listData[["is_penultimate"]] & fs2@elementMetadata@listData[["res_aenmd"]]@listData[["is_last"]],
                       NMDesc_css = fs2@elementMetadata@listData[["res_aenmd"]]@listData[["is_cssProximal"]],
                       NMDesc_long = fs2@elementMetadata@listData[["res_aenmd"]]@listData[["is_407plus"]],
                       key = fs2@elementMetadata@listData[["id"]]
              )

#fs_NMD_result
sum(aenmd_NMD$NMDesc_long == fs_NMD_result$NMDesc_long)
sum(aenmd_NMD$NMDesc_long == fs_NMD_result$NMDesc_long) / 77169
sum(aenmd_NMD$NMDesc_can == fs_NMD_result$NMDesc_can) / 77169
sum(aenmd_NMD$NMDesc_css == fs_NMD_result$NMDesc_css) / 77169