# Wrap the pipeline

library(rentrez)
library(Biostrings)
library(rtracklayer)
library(httr)
library(foreach)
library(doParallel)
library(changepoint)
library(tidyverse)
library(patchwork)
library(cowplot)
library(ggpubr)


GCPeak <- function(accession_number, windows) {
  parts <- unlist(strsplit(accession_number, ""))
  first_part <- paste(parts[1:3], collapse = "")
  second_part <- paste(parts[5:length(parts)], collapse = "")
  fasta_url <- paste0("https://ftp.ncbi.nlm.nih.gov/genomes/all/",
                      first_part,"/",
                      substr(second_part, 1, 3), "/",
                      substr(second_part, 4, 6), "/",
                      substr(second_part, 7, 9), "/")
  # 访问FTP目录
  res <- GET(fasta_url)
  # 解析FTP目录内容
  content <- content(res, "text")
  lines <- strsplit(content, "\n")[[1]]
  # 找到唯一的下一层链接
  next_layer_links <- grep("href=\"", lines, value = TRUE)
  next_layer_links <- gsub(".*href=\"([^\"]+)\".*", "\\1", next_layer_links)
  next_layer_links <- next_layer_links[!grepl("\\.\\.\\/|^\\/", next_layer_links)]
  next_layer_url <- paste0(fasta_url,next_layer_links[[1]])
  # 访问下一层链接
  res <- GET(next_layer_url[[1]])
  next_layer_url[[1]]
  filename <- gsub("\\/","_genomic.fna.gz",next_layer_links[[1]])
  GCAname <- gsub("\\/","",next_layer_links[[1]])
  # 构建完整的下载链接
  genomic_fna_url <- paste0(next_layer_url[[1]], filename)
  if (!dir.exists("GCPeak")) {
    dir.create("GCPeak")
  }
  if (!dir.exists("GCPeak/download")) {
    dir.create("GCPeak/download")
  }
  download.file(genomic_fna_url, destfile = paste0("GCPeak/download/",filename),method = "auto")
  # 定义下载函数
  # download_with_resume <- function(url, destfile) {
  #   # 打开文件准备写入
  #   h <- new_handle()
  #   handle_setopt(h, url = url, resume_from_large = file.info(destfile)$size)
  #   curl_download(url, destfile, handle = h)
  # }
  # download_with_resume(genomic_fna_url, paste0("GCPeak/download/",filename))
  chr2acc_url <- paste0(next_layer_url[[1]],GCAname,"_assembly_structure/Primary_Assembly/assembled_chromosomes/chr2acc")
  filename_chr2acc <- paste0(GCAname,".","chr2acc")
  download.file(chr2acc_url, destfile = paste0("GCPeak/download/",filename_chr2acc))

  # 读取基因组序列
  genome <- readDNAStringSet(paste0("GCPeak/download/",filename))
  df_chr2acc <- read.csv(paste0("GCPeak/download/",filename_chr2acc),sep = "\t")
  chr_acc <- df_chr2acc$Accession.version
  names(genome) <- gsub(" .*","",names(genome))
  genome <- genome[names(genome)%in%chr_acc]
  # 计算GC含量的函数
  calculate_gc_content <- function(seq, window_size) {
    num_seq <- as.numeric(as.vector(seq))
    gc_content <- rollapply(num_seq, width = window_size, by = window_size,
                            FUN = function(x) sum(x %in% c(2, 3)) / window_size)
    return(gc_content)
  }

  # num_cores <- detectCores(8)   # 使用可用核心数减一的核心数
  cl <- makeCluster(8)
  registerDoParallel(cl)
  windows = 20000
  # 用100,000滑窗计算每条染色体的GC含量，并生成结果dataframe
  gc_results <- foreach(chr = names(genome), .combine = rbind, .packages = c("Biostrings", "zoo")) %dopar% {
    seq <- genome[[chr]]
    gc_content <- calculate_gc_content(seq, windows)
    curdf <- data.frame(
      chrom = chr,
      start = seq(1, length(seq) - windows + 1, by = windows),
      GC = gc_content)
    return(curdf)
  }
  stopCluster(cl)

  # gc_results |> ggplot(aes(Start,GC_Content)) +geom_point() + facet_wrap(~Chromosome,scales = "free")
  df1 <- gc_results |>
    select(chrom, start, GC) |>
    filter(GC > 0) |>
    filter(GC < 1) |>
    group_by(chrom) |>
    filter(n() >= 50)
  df_chr2acc$X.Chromosome
  df2 <- df1 |>
    left_join(df_chr2acc |> rename(chrom = Accession.version),
              by = "chrom"
    ) |>
    as.data.frame() |>
    select(-chrom) |>
    rename(chrom = X.Chromosome)
  # custom_order <- c(1:100, "X", "Z", "Y", "W")
  custom_order <- df_chr2acc$X.Chromosome
  df2 <- df2[order(match(df2$chrom, custom_order)), ]
  chrom_list <- df2 |>
    select(chrom) |>
    distinct() |>
    pull(chrom) |>
    as.character()
  grouped_data <- split(df2$GC, df2$chrom)
  result_list <- lapply(grouped_data, function(x) {
    tryCatch(cpt.meanvar(x, method = "BinSeg", penalty = "None", Q = 4, minseglen = 8),
             error = function(e) {
               return(NULL)
             }
    )
  })
  removal <- names(result_list[sapply(result_list, is.null)])
  result_list <- result_list[!sapply(result_list, is.null)]
  chrom_list <- chrom_list[!(chrom_list %in% removal)]

  ###################################### Plot changepoint
  {
    plot_list_cpt <- list()
    for (i in names(result_list)) {
      cur_df <- data.frame(
        x = (1:length(result_list[[i]]@data.set)) * 2 * windows / 1000000,
        y = result_list[[i]]@data.set,
        cpt1 = (result_list[[i]]@cpts[[1]] * 2 * windows) / 1000000,
        cpt2 = (result_list[[i]]@cpts[[2]] * 2 * windows) / 1000000,
        cpt3 = (result_list[[i]]@cpts[[3]] * 2 * windows) / 1000000,
        cpt4 = (result_list[[i]]@cpts[[4]] * 2 * windows) / 1000000,
        cpt5 = (result_list[[i]]@cpts[[5]] * 2 * windows) / 1000000,
        mean1 = result_list[[i]]@param.est$mean[[1]],
        mean2 = result_list[[i]]@param.est$mean[[2]],
        mean3 = result_list[[i]]@param.est$mean[[3]],
        mean4 = result_list[[i]]@param.est$mean[[4]],
        mean5 = result_list[[i]]@param.est$mean[[5]],
        chr = i
      )
      plot_list_cpt[[i]] <- cur_df
    }
    # 合并所有数据框
    combined_df <- bind_rows(plot_list_cpt)
    combined_df$chr <- factor(combined_df$chr, levels = df_chr2acc$X.Chromosome)
    ggplot(combined_df, aes(x, y)) +
      geom_point(size = 2, color = "#666666", fill = "grey", alpha = 1, shape = 21) +
      geom_segment(aes(x = 0, xend = cpt1, y = mean1, yend = mean1), color = "#b71c1c", linewidth = .5) +
      geom_segment(aes(x = cpt1, xend = cpt2, y = mean2, yend = mean2), color = "#b71c1c", linewidth = .5) +
      geom_segment(aes(x = cpt2, xend = cpt3, y = mean3, yend = mean3), color = "#b71c1c", linewidth = .5) +
      geom_segment(aes(x = cpt3, xend = cpt4, y = mean4, yend = mean4), color = "#b71c1c", linewidth = .5) +
      geom_segment(aes(x = cpt4, xend = cpt5, y = mean5, yend = mean5), color = "#b71c1c", linewidth = .5) +
      scale_y_continuous(limits = c(min(df2$GC), max(df2$GC))) +
      theme_classic() +
      theme(legend.position = "none") +
      labs(x = "Chromosome coordinate (Mbp)", y = "GC") +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 12)
      ) +
      facet_wrap(~chr, scales = "free_x")
    ggsave(paste0("GCPeak/cpt.pdf"), plot_wrap_cpt, wi = 12, he = 9, dpi = 150)
  }

  ###################################### Ratio & CV
  # distiance from the tips on each side
  result_df <- data.frame()
  for (i in chrom_list) {
    # tryCatch(
    df_group <- result_list[[i]]
    mean <- df_group@param.est$mean
    cpts <- df_group@cpts

    dat <- tryCatch(
      dat <- data.frame(cpts, mean),
      error = function(e) {
        return(NULL)
      }
    )
    if (!is.null(dat)) {
      seg.1 <- dat[1, ]$cpts
      seg.2 <- dat[2, ]$cpts
      seg.4 <- dat[4, ]$cpts
      # seg.5 is whole length of the chr
      seg.5 <- dat[5, ]$cpts
      seg.1.ratio <- seg.1 / seg.5
      seg.5.ratio <- (seg.5 - seg.4) / seg.5
      seg.1.mean <- dat[1, ]$mean
      seg.2.mean <- dat[2, ]$mean
      seg.4.mean <- dat[4, ]$mean
      seg.5.mean <- dat[5, ]$mean

      # Get chr
      # length.chr <- (length(result_list[[i]]@data.set) * 100000) + 20000
      cur_len <- chr_len$chr_len[chr_len$chrom==i]

      if (seg.1.mean > seg.2.mean && seg.1.ratio < .15) {
        left.ratio <- seg.1.ratio
      } else {
        left.ratio <- NA
      }
      if (seg.5.mean > seg.4.mean && seg.5.ratio < .15) {
        right.ratio <- seg.5.ratio
      } else {
        right.ratio <- NA
      }
      result_df <- rbind(result_df,
                         data.frame(
                           Chr = i,
                           left.ratio = left.ratio, right.ratio = right.ratio,
                           chr_len=cur_len))
    } else {
      print("pass")
    }
  } # Each chromosome

  # Count the ratio and sd
  cur_df <- result_df
  # mean and sd of left or right ratio in each species
  result_df_mean_1 <- cur_df |>
    # select(file_name, left.ratio) |>
    filter(!is.na(left.ratio)) |>
    # group_by(file_name) |>
    summarise(
      left.ratio.mean = mean(left.ratio)
      # left.ratio.sd = sd(left.ratio)
    )
  result_df_mean_3 <- cur_df |>
    # select(file_name, right.ratio) |>
    filter(!is.na(right.ratio)) |>
    # group_by(file_name) |>
    summarise(
      right.ratio.mean = mean(right.ratio)
      # right.ratio.sd = sd(right.ratio)
    )

  seg_1_mean <- result_df_mean_1$left.ratio.mean
  seg_3_mean <- result_df_mean_3$right.ratio.mean

  result_mean_df_list <- data.frame()
  for (x in chrom_list) {
    max_row <- df2[df2$chrom == x, ]
    chr.max.length <- max(max_row$start)

    # Get the max GC value and its position (start)
    # Rationale: chose the highest GC without outliers
    {
      if (!is_empty(target_row_1) && !is_empty(target_row_3)) {
        df.left <- max_row |>
          filter(start < seg_1_mean * chr.max.length)
        df.right <- max_row |>
          filter(start > (1 - seg_3_mean) * chr.max.length)
        # Pearson correlation test
        tryCatch(
          cor.rst.left <- cor.test(df.left$start, df.left$GC, method = "pearson")
          ,error = function(e) {return(NULL)})
        tryCatch(
          cor.rst.right <- cor.test(df.right$start, df.right$GC, method = "pearson")
          ,error = function(e) {return(NULL)})
        # Pearson's coefficient
        r_left <- cor.rst.left$estimate[[1]]
        r_right <- cor.rst.right$estimate[[1]]
        r_left_pval <- cor.rst.left$p.value
        r_right_pval <- cor.rst.right$p.value
      }
      if (is_empty(target_row_1) && !is_empty(target_row_3)) {
        df.left <- max_row |>
          filter(start < seg_3_mean * chr.max.length)
        df.right <- max_row |>
          filter(start > (1 - seg_3_mean) * chr.max.length)
        # Pearson correlation test
        tryCatch(
          cor.rst.left <- cor.test(df.left$start, df.left$GC, method = "pearson")
          ,error = function(e) {return(NULL)})
        tryCatch(
          cor.rst.right <- cor.test(df.right$start, df.right$GC, method = "pearson")
          ,error = function(e) {return(NULL)})
        # Pearson's coefficient
        r_left <- cor.rst.left$estimate[[1]]
        r_right <- cor.rst.right$estimate[[1]]
        r_left_pval <- cor.rst.left$p.value
        r_right_pval <- cor.rst.right$p.value
      }
      if (!is_empty(target_row_1) && is_empty(target_row_3)) {
        df.left <- max_row |>
          filter(start < seg_1_mean * chr.max.length)
        df.right <- max_row |>
          filter(start > (1 - seg_1_mean) * chr.max.length)
        # Pearson correlation test
        tryCatch(
          cor.rst.left <- cor.test(df.left$start, df.left$GC, method = "pearson")
          ,error = function(e) {return(NULL)})
        tryCatch(
          cor.rst.right <- cor.test(df.right$start, df.right$GC, method = "pearson")
          ,error = function(e) {return(NULL)})
        # Pearson's coefficient
        r_left <- cor.rst.left$estimate[[1]]
        r_right <- cor.rst.right$estimate[[1]]
        r_left_pval <- cor.rst.left$p.value
        r_right_pval <- cor.rst.right$p.value
      }
      if (is_empty(target_row_1) && is_empty(target_row_3)) {
        df.left.mean <- -1
        df.right.mean <- -1
        # Pearson correlation test
        # tryCatch(
        #   cor.rst.left <- cor.test(df.left$start, df.left$GC, method = "pearson")
        #   ,error = function(e) {return(NULL)})
        # tryCatch(
        # cor.rst.right <- cor.test(df.right$start, df.right$GC, method = "pearson")
        #   ,error = function(e) {return(NULL)})
        # Pearson's coefficient
        r_left <- -99
        r_right <- -99
        r_left_pval <- -99
        r_right_pval <- -99
      } # no cpt, then discard
    }
    # Get the CV of GC conetent of whole chr (Remove outlier)
    GC.mean <- max_row$GC |> mean()
    GC.sd <- max_row$GC |> sd()
    GC.cv <- GC.sd / GC.mean
    result_mean_df_list <- rbind(
      result_mean_df_list,
      data.frame(
        # file_name = i,
        chr = x,
        r_left = r_left,
        # r_left_pval = r_left_pval,
        r_right = r_right,
        # r_right_pval = r_right_pval,
        GC.mean=GC.mean,
        # GC.sd=GC.sd,
        GC.cv=GC.cv,
        chr_len = chr.max.length
        # species = cur_species
      )
    )

  } # Each chromosome

  ###################################### Plot final
  result_pearson <- result_mean_df_list |>   mutate(left = ifelse( r_left < (-0.2) , 1, 0),
                                                    right = ifelse(r_right > 0.2, 2, 0)) |>
    mutate(spike = right + left)|>
    mutate(
      quantification = ifelse(
        spike == 3, {
          max_abs <- pmax(abs(r_left), abs(r_right))
          (max_abs)
        },
        ifelse(spike == 2, abs(r_right),
               ifelse(spike == 1, abs(r_left), 0))
      )
    )
  dat_text <- data.frame(
    label = paste0("Ratio:",round(result_pearson$quantification,2)," CV:",round(result_pearson$GC.cv,2),
                   " GC:",round(result_pearson$GC.mean,2)),
    chrom = df2$chrom |> unique())

  df2 |> ggplot(aes(start/1000000,GC)) +
    geom_point(size = 1.5, color = "#666666", fill = "grey", alpha = 1, shape = 21) +
    facet_wrap(~chrom, scales = "free_x") +
    geom_text(data = dat_text,
              mapping = aes(x = 0, y = max(df2$GC), label = label),
              hjust   = 0,vjust   = 0) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(x = "Chromosome coordinate (Mbp)", y = "GC",
         title = paste0("Overall Ratio:", round(mean(result_pearson$quantification),2)," CV:",round(mean(result_pearson$GC.cv),2), " GC:",round(mean(result_pearson$GC.mean),2))
    ) +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 12, color = "black"),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 12),
      strip.background = element_blank(),
      strip.text = element_text(hjust = 0, size = 12)
    )  +
    scale_y_continuous(expand = expansion(mult = c(.01,.1)))
  ggsave("GCPeak/cpt.pdf",wi=12,he=9)
}

GCPeak("GCA_037355615.1",1000000)

