# fourier
println("start compile...")
using PyPlot
pygui(true)
gcf()
using FFTW
using Printf
#∭∭∭∭∭∭∭∭∭ 常數  ∭∭∭∭∭∭∭∭∭∭
const τ = 8
const time_length = 500
const down_range = -6π
const up_range = 6π
const f(x) = exp(-x^2*τ)
#∭∭∭∭∭∭∭∭∭∭∭∭∭∭∭∭∭∭∭∭∭
function make_shift_freq(n,step)
    #f = [0, 1, ..., (n-1)/2, -(n-1)/2, ..., -1] / (d*n)   if n is odd
    if (n % 2 ==1)
    out = collect(Iterators.flatten([collect(0:(n-1/2)),collect(-(n-1)/2:-1)])) / (step * n)
    elseif (n % 2 ==0)
    #f = [0, 1, ...,   n/2-1,     -n/2, ..., -1] / (d*n)   if n is even
    out = collect(Iterators.flatten([collect(0:n/2), collect(-n/2:-1)]))/(step * n)
    end
    #shift freq time
    return fftshift(out)
end
function find_FWHM(args,time,ff = false)
    # after fourier
    if (ff)
        suppose =  √(2*log(2))/π*√τ
        @printf("理論值: %9.7f   傅立葉後",suppose)
        f = make_shift_freq(length(time),abs(time[1] - time[2]))
    # before fourier
    else
        suppose =  √(2*log(2))/√τ
        @printf("理論值: %9.7f   傅立葉前",(suppose))
        f = time
    end
    # Intensity
    IF = map(x->x^2,args)
    # find the half max 半高
    get = 0.5 * maximum(IF)
    # 將矩陣的每個值都剪掉剛剛找到的最大值的一半
    IF = map(x->get - x,IF)
    # 做絕對值
    IF = map(abs,IF)
    # 找到離 0 對小那點的x 取絕對值值 *2 就是 fwhm (因為高斯是對稱的)
    # 如果不是高斯要換算法
    fwhm = 2 * abs( f[ findmin(IF)[2] ] )
    @printf("半高寬:   %9.7f 誤差: %7.5f %s",fwhm, (fwhm - suppose)/suppose*100,"%")
    println()
    return fwhm
end
function start()
    # 時間矩陣
    time = collect(range(down_range,up_range,length = time_length))
    println("時間矩陣大小:   ", length(time))
    # 製造 高斯的函數
    y = map(f,time)
    # 去找半高寬
    yFWHM = find_FWHM(y,time)
    # 做 傅立葉 和 shift 後去做 絕對值
    y_fft = map(abs,fftshift(fft(y)))
    y_fft  = y_fft / findmax(y_fft)[1]
    # 去找 傅立葉後的 半高寬
    y_fftFWHM = find_FWHM(y_fft,time,true)
    # 相乘
    out = y_fftFWHM * yFWHM
    # 看誤差
    @printf("理論值: 0.441       FWHM :  %9.7f 誤差: %7.5f %s",out, (out - 0.441)/0.441*100,"%")
    # 畫圖
    subplot(211)
    title("before fourier")
    plot(time,y)
    subplot(212)
    title("after fourier")
    fre = make_shift_freq(length(time),abs(time[1] - time[2]))
    #  fre[1:end-1]
    plot(fre[1:end-1],y_fft,"r")
    println()
    print("所花時間:")
end
@time start()
